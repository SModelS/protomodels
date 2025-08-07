#!/usr/bin/env python3

""" module that contains stuff on top of SModelS """

def suppressNoValidDecay():
    """ we may currently have particles with no valid decay.
    they are set to stable. FIXME: should be explicitly set to stable. """
    from smodels.base.smodelsLogging import logger, logging

    class SuppressNoValidDecay(logging.Filter):
        def filter(self, record):
            return "No valid decay found for " not in record.getMessage()

    logger.addFilter(SuppressNoValidDecay())
