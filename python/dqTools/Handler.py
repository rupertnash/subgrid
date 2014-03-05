"""Adds the ability to Controllers to checkpoint upon receipt of a
signal, by default SIGTERM. There can be only one instance of Handler.

Your class must subclass from both Handler.Handler and
Controller.Controller, furthermore Handler must be the leftmost of the
superclasses. E.g.,

class MyController(Handler.Handler, Controller.Controller):
    pass

"""

import sys
import signal
import weakref

class Handler(object):
    # The commented out bits were from an attempt to relax the
    # singleton requirement.
    instance = None#weakref.WeakValueDictionary()
    #i = 0
    def __new__(cls, *args, **kwargs):
        if Handler.instance is None or Handler.instance() is None:
            self = object.__new__(cls, *args, **kwargs)
            self._signalled = False
            Handler.instance = weakref.ref(self)
        else:
            raise RuntimeError('This is a singleton class that inherits from Handler!')

        #cls.instances[cls.i] = self
        #cls.i += 1
        return self

    @classmethod
    def handle(cls, signum, frame):
        #for con in cls.instances.values():
        if cls.instance is not None:
            con = cls.instance()
            if con is not None:
                # easiest would be to call con.checkpoint() but we can't
                # just stop part way through a timestep, need to let this
                # one finish, THEN checkpoint (& die). Instead, set
                # con._signalled to True and alter isCPStep to deal.
                con._signalled = True
                con.log('Received signal to checkpoint and die')
        return

    def isCheckpointStep(self):
        if self._signalled:
            return True

        return super(Handler, self).isCheckpointStep()

    def checkpoint(self):
        signalled = self._signalled
        # We really don't want to be saving the state with this flag set:
        # as soon as it runs on resume it will checkpoint & die!
        self._signalled = False
        super(Handler, self).checkpoint()
        if signalled:
            self.log('Terminating due to signal... bye!')
            sys.exit(0)
            pass

        return

    pass

signal.signal(signal.SIGTERM, Handler.handle)


