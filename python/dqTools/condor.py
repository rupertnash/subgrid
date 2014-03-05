"""Guess if we're running as a condor job.

This is a hack-- test carefully on your system!
"""
import os

if len([k for k in os.environ if k.startswith('_CONDOR_ANCESTOR')]):
    onCondor = True
else:
    onCondor = False

