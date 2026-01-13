#!/usr/bin/env python3

""" Class for all things around logging """

import time, os
from typing import Union
from ptools import helpers
from smodels_utils.helper.terminalcolors import *

__all__ = [ "LoggerBase" ]

class LoggerBase:
    __slots__ = [ "walkerid", "module", "logdir" ]

    def __init__ ( self, walkerid : Union[str,int] = 0 ):
        """ instantiate the logger class with a walkerid """
        self.walkerid = walkerid
        self.countLogs = {}
        self.printLogMessages = False
        self.logdir = "logs/"
        # self.printHigherThan = "critical"
        module = str(type(self)).replace("<class '","").replace("'>","")
        p1 = module.find(".")
        if module.count(".")==2:
            p2 = module.rfind(".")
            self.module = module[p1+1:p2]
        else:
            self.module = module[p1+1:]
        helpers.mkdir ( self.logdir )

    def logThrice ( self, *args ):
        """ a method for repetitive log msgs. issue them only three times
        """
        txt = " ".join(map(str,args))
        if not txt in self.countLogs:
            self.countLogs[txt]=0
        if self.countLogs[txt]<3:
            self.log (  *args )
        if self.countLogs[txt]==3:
            self.log ( "(quenching repeating log messages)" )
        self.countLogs[txt]+=1

    def error ( self, *args ):
        self.highlight ( "error", *args )

    def warn ( self, *args ):
        self.highlight ( "warn", *args )
        
    def warning ( self, *args ):
        self.highlight ( "warn", *args )

    def info ( self, *args ):
        """ logging to file, but also write to screen """
        self.log ( *args )
        if self.verbose > 1:
            print ( f"[expResModifier] {' '.join(map(str, args))}" )

    def highlight ( self, msgType : str = "info", *args ):
        """ logging, hilit """
        col = GREEN
        if msgType.lower() in [ "error", "red" ]:
            col = RED
        elif msgType.lower() in [ "warn", "warning", "yellow" ]:
            col = YELLOW
        elif msgType.lower() in [ "green", "info" ]:
            col = GREEN
        else:
            self.highlight ( "red", "I think we called highlight without msg type" )
        print ( f'{col}[{self.module}:{time.strftime("%H:%M:%S")}] {" ".join(map(str,args))}{RESET}' )
        self.log ( *args )

    def debug ( self, *args ):
        tmp = list ( args )
        for i,arg in enumerate ( tmp ):
            if type ( arg ) == str:
                tmp[i] = 'DEBUG: ' + arg
                break 
        args = tuple ( tmp )
        self.log ( *args )

    def pprint ( self, *args ):
        """ logging """
        line = ' '.join(map(str,args))
        if not line in self.countLogs:
            self.countLogs[line]=0
        self.countLogs[line] += 1
        if self.countLogs[line]==3:
            self.pprint ( f"skipping repeating messages" )
            return
        if self.countLogs[line]>3:
            return
        print ( f"[{self.module}:{self.walkerid}] {line}" )
        self.prevMessage = line
        self.log ( *args )

    def cprint ( self, color, *args ):
        """ logging, colored version """
        line = ' '.join(map(str,args))
        from smodels_utils.helper.terminalcolors import colordict, RESET
        if not line in self.countLogs:
            self.countLogs[line]=0
        self.countLogs[line] += 1
        if self.countLogs[line]==3:
            self.pprint ( f"skipping repeating messages" )
            return
        if self.countLogs[line]>3:
            return
        print ( f"[{self.module}:{self.walkerid}] {colordict[color]}{line}{RESET}" )
        self.prevMessage = line
        self.log ( *args )

    def log ( self, *args ):
        """ logging to file """
        helpers.mkdir ( self.logdir )
        ctr = 0
        while True:
            try:
                with open( f"{self.logdir}/walker_{self.walkerid}.log", "a" ) as f:
                    f.write ( f'[{self.module}-{time.strftime("%H:%M:%S")}] {" ".join(map(str,args))}\n' )
                if self.printLogMessages:
                    print ( f'[{self.module}-log] {" ".join(map(str,args))}' )

                return
            except OSError as e:
                # lets try a few times, we are using network file systems,
                # the network might be acting out
                ctr+=1
                time.sleep ( ctr**2 )
                if ctr > 10:
                    raise e

if __name__ == "__main__":
    logger = LoggerBase ( 0 )
    logger.pprint ( "what now!" )
