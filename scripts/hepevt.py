#! /bin/env python
#=================================================================================
#
# Name: hepevt.py
#
# Purpose: CGI script for serving HepEvt or HepMC events.
#
# Created: 11-May-2022  H. Greenlee
#
# Usage:
#
# hepevt.py <arguments>
#
# CLI arguments (for testing):
#
# -h|--help               - Print help message.
# -f|--file <hepevt-file> - HepEvt or HepMC file.
# --format <format>       - File format, "hepevt" or "hepmc" (default "hepevt").
# -s|--stream <stream>    - Event counter stream name (default='default').
# -e|--event <event>      - Specify event number.
# -r|--reset              - Reset event counter for the specified stream.
# --sleep <secs>          - Sleep time per event (default 0)
# --min_event <event>     - Minimum event number.
# --max_event <event>     - Maximum event number.
#
# CGI arguments:
#
# file      - Name of hepevt or hepmc file.
# format    - File format, "hepevt" or "hepmc" (default "hepevt").
# stream    - Event counter stream name (default='default').
# event     - Specify event number.
# reset     - Event counter reset flag (reset if "1").
# sleep     - Sleep time per event.
# min_event - Minimum event number.
# max_event - Maximum event number.
#
# Usage notes.
#
# 1.  This script is python 2/3 agnostic.  When invoked on the web server, it uses
#     the system python (current python 2.7 on the MicroBooNE web server).  It also
#     works using larsoft environment pythons, whether python 2 or 3.
#
# 2.  This script accepts command line arguments or CGI arguments embedded in a url.
#
# 3.  The HepEvt or HepMC file name should be specified as a plain file or relative
#     path name (absolute path names are not allowed for security reasons).  If a
#     file specified using a relative path is not found relative to the current
#     directory, this script searches for the file using a search path that consists
#     of the following directories.
#
#     a) $DOCUMENT_ROOT/../data/hepevt
#     b) /web/sites/<letter>/<group>-exp.fnal.gov/data/hepevt   (linux group name)
#     c) /web/sites/<letter>/${GROUP}-exp.fnal.gov/data/hepevt  (group env variable)
#
# 4.  HepMC format can be auto-sensed, making the format argument optional for either
#     HepEvt or HepMC format files.
#
# 5.  The event counter file is located in the same directory as the HepEvt file.
#     This means that this script requires write access to this directory.
#
# 6.  The event counter file contains the event number and byte offset of the most
#     recently read HepEvt event.  If no events have been read, the event counter file
#     doesn't exist.
#
# 7.  The name of event counter file is as follows:
#     <hepevt-file>.<user>.<stream>.txt
#     The user name is only included if using SSO authentication.
#     The stream name is as specified by the CGI or CLI argument.
#
# 8.  If a minimum or maximum event number is specified, returned events are limited
#     to the range min_event <= e < max_event.
#
# 9.  This script uses posix file locking to prevent multiple processes from
#     updating the event counter at the same time.  If the script can't acquire a
#     lock, this script reutrns error 503 (service unavailable) so as to not block
#     the entire server.  In such cases, the client should retry the operation.
#     Resetting the event counter does not make use of file locking.
#
# 10.  The sleep time per event option can be used to reduce rate at which events
#      can be delivered.  It is intended for testing, and is not useful for production.
#
# 11.  If a specific event number is specified, the event counter file is not
#      read or updated.
#
# 12.  This script may return the following error codes.
#      404 - Not found.
#            This error is returned for various errors, including:
#            a) No HepEvent file was specified.
#            b) Specified HepEvent file does not exist.
#            c) Invalud event number or attempt to read past end of file.
#      503 - Service unavailable.
#            Couldn't acquire exclusive lock on event counter.  Client should retry.
#
#
#=================================================================================

from __future__ import print_function
import sys, os, fcntl, time, grp
import cgi
import cgitb
cgitb.enable()


# Global variables

file_format = ''


# Print help.

def help():

    filename = sys.argv[0]
    file = open(filename)

    doprint=0
    
    for line in file.readlines():
        if line.startswith('# hepevt.py'):
            doprint = 1
        elif line.startswith('# Usage notes'):
            doprint = 0
        if doprint:
            if len(line) > 2:
                print(line[2:], end='')
            else:
                print()


# Seek from the current file position to the specified, or greater, event number.
#
# Returns the actual event number, or -1 if no matching event is found.

def seek_event_number(fhepevt, evnum):

    global file_format

    result = -1

    # Read lines until we find a compatible event number.

    while True:
        pos = fhepevt.tell()
        line = fhepevt.readline()
        if line == '':
            result = -1
            break
        words = line.split()

        # Check event number.

        if file_format == 'hepevt':

            # Hepevt format.
            # New events are indicated by a line with two words:
            # <evnum> <nparticles>

            if len(words) == 2 and words[0].isdigit() and words[1].isdigit():
                n = int(words[0])
                if n >= evnum:
                    fhepevt.seek(pos)
                    result = n
                    break

            # If this looks like a hepmc event header, switch the format.

            elif len(words) == 4 and words[0] == 'E' and words[1].isdigit() and \
               words[2].isdigit() and words[3].isdigit():
                file_format = 'hepmc'
                fhepevt.seek(pos)
                continue


        elif file_format == 'hepmc':

            # Hepmc format.
            # New events are indicated by a line with four words:
            # E <evnum> <nvertices> <nparticles>

            if len(words) == 4 and words[0] == 'E' and words[1].isdigit() and \
               words[2].isdigit() and words[3].isdigit():
                n = int(words[1])
                if n >= evnum:
                    fhepevt.seek(pos)
                    result = n
                    break

        else:

            # Unknown format.

            break


    # Done.

    return result


# Function to seek to the next event, as specified in the event counter file.
#
# If successful, this function will increment and update the event counter.
#
# Returns an error code, which can be any of the following.
#
# 0   - Success.
# 404 - Couldn't find event.
# 503 - Couldn't obtain an exclusive lock on event counter file.

def seek_next_event(fhepevt, event_counter_path, min_event, max_event, sleep_time):

    # Open the event counter file for updating.
    # This open mode will create the file if it doesn't exist, but won't truncate it if 
    # it does exist.

    flock = open(event_counter_path, 'a+')
    flock.seek(0, 2)           # Seek to end of file (necessary in python 2.7)

    # Attempt to obtain an exclusive lock.
    # Return error if unsuccessful.

    lockok = False
    try:
        fcntl.lockf(flock, fcntl.LOCK_EX | fcntl.LOCK_NB)  # Exclusive non-blocking.
        lockok = True
    except:
        lockok = False
    if not lockok:
        flock.close()
        return 503

    # Read existing event number, if any.
    # Seek to last read event position.

    evnum = -1
    if flock.tell() > 0:
        flock.seek(0)
        line = flock.readline()
        words = line.split()
        evnum = int(words[0])
        pos = int(words[1])
        fhepevt.seek(pos)
    else:
        evnum = -1

    # Seek to the next event number.

    seek_event = evnum + 1
    if seek_event < min_event:
        seek_event = min_event
    new_evnum = seek_event_number(fhepevt, seek_event)
    if new_evnum < 0 or (max_event > 0 and new_evnum >= max_event):
        flock.close()
        return 400

    # Store the newly found event number and position back in the event counter file.

    flock.seek(0)
    flock.truncate()
    flock.write('%d %d\n' % (new_evnum, fhepevt.tell()))

    # Close event counter file.  This will also release the lock.

    if sleep_time > 0.:
        time.sleep(sleep_time)
    flock.close()

    # Done

    return 0


# Read and print out one event from file.

def read_event(f):

    global file_format

    if file_format == 'hepevt':

        # Hepevt format.
        # Read and print event header line.

        line = f.readline()
        print(line, end='')
        words = line.split()
        npart = int(words[1])

        # Read and print particles.

        while npart > 0:
            line = f.readline()
            print(line, end='')
            npart -= 1

    elif file_format == 'hepmc':

        # Hepmc format.
        # Read and print event header line.

        line = f.readline()
        print(line, end='')

        # Read and print remaining lines until we get to the next event or footer or end-of-file.

        while True:
            pos = f.tell()
            line = f.readline()
            if line == '':
                break
            words = line.split()
            if words[0] == 'E':
                f.seek(pos)
                break
            if line.startswith('HepMC::'):
                break
            print(line, end='')

    # Done.

    return


# Construct HepEvt search path.

def search_path():

    result = []

    # Based on $DOCUMENT_ROOT
    # This is the main method on the web server.

    if 'DOCUMENT_ROOT' in os.environ:
        ddir = os.environ['DOCUMENT_ROOT']
        pdir = os.path.dirname(ddir)
        dir = os.path.join(pdir, 'data/hepevt')
        if os.path.isdir(dir) and not dir in result:
            result.append(dir)

    # Based on $EXPERIMENT.

    if 'EXPERIMENT' in os.environ:
        exp = os.environ['EXPERIMENT']
        dir = '/web/sites/%s/%s-exp.fnal.gov/data/hepevt' % (exp[0], exp)
        if os.path.isdir(dir) and not dir in result:
            result.append(dir)

    # Based on linux group name.

    exp = grp.getgrgid(os.getgid())[0]
    dir = '/web/sites/%s/%s-exp.fnal.gov/data/hepevt' % (exp[0], exp)
    if os.path.isdir(dir) and not dir in result:
        result.append(dir)

    # Done

    return result


# Main procedure.

def main(argv):

    global file_format

    # Extract arguments.

    hepevt_file_name = ''
    file_format = 'hepevt'
    stream_name = 'default'
    evnum = None
    reset = None
    sleep_time = 0
    min_event = 0
    max_event = 0

    # Parse command line arguments.

    args = argv[1:]
    while len(args) > 0:
        if args[0] == '-h' or args[0] == '--help':

            # Help.
            
            help()
            return 0

        elif args[0] == '-r' or args[0] == '--reset':

            # Reset event counter flag.
            
            reset = '1'
            del args[0]

        elif len(args) > 1 and (args[0] == '-f' or args[0] == '--file'):

            # HepEvt file.
            
            hepevt_file_name = args[1]
            del args[0:2]

        elif len(args) > 1 and args[0] == '--format':

            # File format, hepevt or hepmc.

            file_format = args[1]
            del args[0:2]

        elif len(args) > 1 and (args[0] == '-s' or args[0] == '--stream'):

            # Stream name
            
            stream_name = args[1]
            del args[0:2]

        elif len(args) > 1 and (args[0] == '-e' or args[0] == '--event'):

            # Event number
            
            evnum = int(args[1])
            del args[0:2]

        elif len(args) > 1 and args[0] == '--sleep':

            # Sleep time per event
            
            sleep_time = float(args[1])
            del args[0:2]

        elif len(args) > 1 and args[0] == '--min_event':

            # Minimum event number
            
            min_event = int(args[1])
            del args[0:2]

        elif len(args) > 1 and args[0] == '--max_event':

            # Maximum event number
            
            max_event = int(args[1])
            del args[0:2]

        elif args[0][0] == '-':

            # Unknown option.

            print('Unknown option %s' % args[0])
            return 1
            
        else:

            # Arguments.
            
            print('Unknown argument %s' % args[0])
            return 1

    # Parse CGI arguments.

    args = cgi.FieldStorage()
    if 'file' in args:
        hepevt_file_name = args['file'].value
    if 'format' in args:
        file_format = args['format'].value
    if 'stream' in args:
        stream_name = args['stream'].value
    if 'event' in args:
        evnum = int(args['event'].value)
    if 'reset' in args:
        reset = args['reset'].value
    if 'sleep' in args:
        sleep_time = float(args['sleep'].value)
    if 'min_event' in args:
        min_event = int(args['min_event'].value)
    if 'max_event' in args:
        max_event = int(args['max_event'].value)

    # It is an error of HepEvt file is not specified.

    if hepevt_file_name == '':
        print('Content-type: text/plain')
        print('Status: 404 Not Found')
        print('')
        print('No HepEvt file specified.')
        return 0

    # Reject absolute file paths.

    if hepevt_file_name[0] == '/' or hepevt_file_name[0] == '.':
        print('Content-type: text/plain')
        print('Status: 404 Not Found')
        print('')
        print('Absolute path not allowed.')
        return 0

    # Locate HepEvt file.
    # Return error 404 if not found.

    hepevt_file_path = hepevt_file_name
    if not os.path.exists(hepevt_file_path):
        for dir in search_path():
            hepevt_file_path = os.path.join(dir, hepevt_file_name)
            if os.path.exists(hepevt_file_path):
                break
    if not os.path.exists(hepevt_file_path):
        print('Content-type: text/plain')
        print('Status: 404 Not Found')
        print('')
        print('HepEvt file %s does not exist.' % hepevt_file_name)
        return 0

    # Check file format.

    if file_format != 'hepevt' and file_format != 'hepmc':
        print('Content-type: text/plain')
        print('Status: 404 Not Found')
        print('')
        print('Invalid file format %s.' % file_format)
        return 0

    # Construct path of the event counter file.
    # At this point, this file may or may not exist.

    event_counter_path = ''
    if 'SSO_USERID' in os.environ:
        event_counter_path = '%s.%s.%s.next' % (hepevt_file_path, os.environ['SSO_USERID'], stream_name)
    else:
        event_counter_path = '%s.%s.next' % (hepevt_file_path, stream_name)

    # Reset event counter?

    if reset == '1':
        if evnum == None:
            if os.path.exists(event_counter_path):
                os.remove(event_counter_path)
        else:
            f = open(event_counter_path, 'w')
            f.write('%d 0\n' % evnum)
            f.close()
        print('Content-type: text/plain')
        print('Status: 200 OK')
        print('')
        if evnum == None:
            print('Reset event counter.')
        else:
            print('Reset event counter to %d.' % evnum)
        return 0

    # Open HepEvt file.

    fhepevt = open(hepevt_file_path)

    # Seek to the desired event, either by event number or event counter.

    status = 0
    if evnum == None:

        # Seek to next event.

        status = seek_next_event(fhepevt, event_counter_path, min_event, max_event, sleep_time)

    else:

        # Seek to event number

        new_evnum = seek_event_number(fhepevt, evnum)
        if new_evnum < 0:
            status = 400

    # Check for errors.

    if status == 400:
        print('Content-type: text/plain')
        print('Status: 400 Seek Error')
        print('')
        print('Event not found.')
        return 0
    elif status == 404:
        print('Content-type: text/plain')
        print('Status: 404 Not Found')
        print('')
        print('Event not found.')
        return 0
    elif status == 503:
        print('Content-type: text/plain')
        print('Status: 503 Service unavailable')
        print('')
        print('Access to event counter temporarily unavailable.')
        return 0
    elif status != 0:
        print('Content-type: text/plain')
        print('Status: %d' % status)
        print('')
        return 0

    # Don't expect errors from here on.
    # Print header.

    print('Content-type: text/plain')
    print('Status: 200 OK')
    print('')

    # Print selected event.

    read_event(fhepevt)

    # Done.

    return 0


# Invoke main program.

if __name__ == "__main__":
    sys.exit(main(sys.argv))
