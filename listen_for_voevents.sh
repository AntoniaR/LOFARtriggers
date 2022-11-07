#!/bin/sh

#source /home/antoniar/virtenv/bin/activate
#source /home/antoniar/TraP_Develop_Env/bin/activate
source /home/antoniar/python3.6_env/bin/activate

# Start an instance of Comet and connect to 4PiSky broker.
# VOEvent packets received are passed via stdin to the $HANDLER script.


# The 'local ivorn' is the string identifying this machine on the VOEvent network.
# It's not very important unless you want to *broadcast* events, but it has
# to conform to a standard format all the same, something like e.g.
# ivo://fpstoolstest/some_random_machine
LOCALIVORN=ivo://org.lofartkp/struis

# This is the script that comet will pass every VOEvent to, for processing:
#HANDLER=./sgrb_swift.py
HANDLER=./LOFAR_triggers.py

# The remote broker that we will subscribe to. The '--remote' flag can be
# repeated with different addresses if you want to  listen to multiple remotes.
REMOTE=voevent.4pisky.org


# Put it all together:
/usr/bin/env twistd -n comet \
    --verbose \
    --receive \
    --local-ivo=$LOCALIVORN \
    --remote=$REMOTE \
    --cmd=$HANDLER



