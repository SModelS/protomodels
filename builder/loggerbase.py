#!/usr/bin/env python3

""" Class for all things around logging """

from colorama import Fore as ansi
import time

__all__ = [ "LoggerBase" ]

class LoggerBase:
    __slots__ = [ "walkerid", "module" ]

    def __init__ ( self, walkerid : int = 0 ):
        """ instantiate the logger class with a walkerid """
        self.walkerid = walkerid
        module = str(type(self)).replace("<class '","").replace("'>","")
        p1 = module.find(".")
        if module.count(".")==2:
            p2 = module.rfind(".")
            self.module = module[p1+1:p2]
        else:
            self.module = module[p1+1:]

    def error ( self, *args ):
        self.highlight ( "error", *args )

    def info ( self, *args ):
        """ logging to file, but also write to screen """
        self.log ( *args )
        if self.verbose > 1:
            print ( "[expResModifier] %s" % ( " ".join(map(str,args)) ) )


    def highlight ( self, msgType : str = "info", *args ):
        """ logging, hilit """
        col = ansi.GREEN
        if msgType.lower() in [ "error", "red" ]:
            col = ansi.RED
        elif msgType.lower() in [ "warn", "warning", "yellow" ]:
            col = ansi.YELLOW
        elif msgType.lower() in [ "green", "info" ]:
            col = ansi.GREEN
        else:
            self.highlight ( "red", "I think we called highlight without msg type" )
        print ( f'{col}[{self.module}:{time.strftime("%H:%M:%S")}] {" ".join(map(str,args))}{ansi.RESET}' )
        self.log ( *args )

    def debug ( self, *args ):
        pass

    def pprint ( self, *args ):
        """ logging """
        print ( f"[{self.module}] {' '.join(map(str,args))}" )
        self.log ( *args )

    def log ( self, *args ):
        """ logging to file """
        with open( f"walker{self.walkerid}.log", "a" ) as f:
            f.write ( f'[{self.module}-{time.strftime("%H:%M:%S")}] {" ".join(map(str,args))}\n' )

if __name__ == "__main__":
    logger = LoggerBase ( 0 )
    logger.pprint ( "what now!" )
