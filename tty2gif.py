#!/usr/bin/env python
"""Script to convert ttyrec file to gifs on OSX with Terminal.app."""

import argparse
import os
import sys
import struct
import glob
import time
import subprocess

FrameHead = struct.Struct('iii')

parser = argparse.ArgumentParser(description='Convert ttyrec data to gifs.')
parser.add_argument('action', choices=['replay', 'inspect', 'output'],
                   help='Action to be performed on the ttyrec binary.')
parser.add_argument('input', type=str, help='ttyrec filename to be processed.')
parser.add_argument('--output', type=str, help='Output gif filename')
parser.add_argument('--factor', type=int, default=1, help='Speedup factor.')
parser.add_argument('--cropframe', type=bool, default=False,
                    help='Whether to crop the window frame or not.')

SKIP_THREASHOLD = 0.005
SKIP_LIMIT = 5

def _take_screenshot(winid, filename):
    """Take screenshot of current terminal."""
    # capture current frame in png
    os.system('screencapture -l %s -o %s' % (winid, filename))

def replay(payload, delay, skip=False):
    """No-op placeholder for replay action.

    The frame rate of replay reflects the delay of ttyrec with `--factor`
    parameter to adjust the speed.
    """
    sys.stdout.write(payload)
    sys.stdout.flush()
    time.sleep(delay)

def output(payload, delay, skip=False):
    """Write current frame to files.

    TODO: Different from replay, the frame rate of capturing is constant to
    preseve the original speed in final output gif.

    """
    output.winid = getattr(output, 'winid', None)
    if not output.winid:
        # NOTE: could only work with "Terminal.app", 'iTerm" doesn't
        # return the correct window id.
        output.winid = int(subprocess.check_output([
            'osascript', '-e', 'tell app "Terminal" to id of window 1']))

    if not skip:
        output.count = getattr(output, 'count', -1) + 1
        png_file = 'step_%03d.png' % output.count
        _take_screenshot(output.winid, png_file)

    sys.stdout.write(payload)
    sys.stdout.flush()
    time.sleep(0.005)

def inspect(payload, delay, skip=False):
    """Action to print out the frame data for inspection."""
    n = len(payload)
    print('%8.4f %4d %s' % (delay, n, payload[:40]))

if __name__ == '__main__':
    args = parser.parse_args()
    try:
        action = globals()[args.action]
        factor = args.factor
    except KeyError:
        print('Action %s is not defined as method.' % args.action)
        exit(1)

    with open(args.input) as script:
        basetime = None
        prevtime = None
        skip = False
        nskipped = 0
        while True:
            data = script.read(FrameHead.size)
            if not data:
                break
            sec, usec, n = FrameHead.unpack(data)
            payload = script.read(n)

            # compute the delay from last frame
            curtime = sec + usec / 1000000.0
            if basetime is None:
                basetime = prevtime = curtime
            delay = (curtime - prevtime) / factor
            prevtime = curtime


            # Skip actions for some
            if (delay <= SKIP_THREASHOLD):
                skip = True
                nskipped += 1
            else:
                skip = False
                nskipped = 0

            if (skip and nskipped > SKIP_LIMIT):
                nskipped = 0
                skip = False

            action(payload, delay, skip)

    # remove transparency and crop image combined

    screenshots = glob.glob('step_*.png')
    if args.cropframe:
        print("cropping off window frame...")
        for filename in screenshots:
            os.system('convert %s -background white -flatten +matte -crop +0+8 -crop +4+0 -crop -4-0 %s' % (filename, filename))

    # print "generating the final gif..."
    os.system('convert -loop 0 %s -layers Optimize %s' % (' '.join(screenshots), 'tty.gif'))
