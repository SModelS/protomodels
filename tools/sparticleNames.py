#!/usr/bin/env python3

"""
.. module:: sparticleNames
        :synopsis: assign sparticle names to pids ( 1000021 <-> ~g or Xg, ... ),
        pids to names, categorizes particles, etc.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
import sys

class SParticleNames:
    """ a class that assigns names to sparticles """

    def initSUSY ( self ):
        """ define the SUSY notation """
        self.susyIDs = { 1000001: "~d_L", 2000001: "~d_R",
            1000002: "~u_L", 2000002: "~u_R",
            1000003: "~s_L", 2000003: "~s_R",
            1000004: "~c_L", 2000004: "~c_R",
            1000005: "~b_1", 2000005: "~b_2",
            1000006: "~t_1", 2000006: "~t_2",
            1000011: "~e_L", 1000012: "~nu_eL",
            1000013: "~mu_L", 1000014: "~nu_muL",
            1000015: "~tau_L", 1000016: "~nu_tauL",
            2000011: "~e_R", 2000012: "~nu_eR",
            2000013: "~mu_R", 2000014: "~nu_muR",
            2000015: "~tau_R", 2000016: "~nu_tauR",
            1000021: "~g", 1000022: "~chi10",
            1000024: "~chi1+", 1000023: "~chi20",
            1000037: "~chi2+", 1000025: "~chi30",
            1000035: "~chi40",
            -1000001: "~d_L^*", -2000001: "~d_R^*",
            -1000002: "~u_L^*", -2000002: "~u_R^*",
            -1000003: "~s_L^*", -2000003: "~s_R^*",
            -1000004: "~c_L^*", -2000004: "~c_R^*",
            -1000005: "~b_1^*", -2000005: "~b_2^*",
            -1000006: "~t_1^*", -2000006: "~t_2^*",
            -1000011: "~e_L^*", -1000012: "~nu_eL",
            -1000013: "~mu_L^*", -1000014: "~nu_muL",
            -1000015: "~tau_L^*",- 1000016: "~nu_tauL",
            -2000011: "~e_R^*", -2000012: "~nu_eR",
            -2000013: "~mu_R^*", -2000014: "~nu_muR",
            -2000015: "~tau_R^*", -2000016: "~nu_tauR",
            -1000021: "~g^*", -1000022: "~chi10",
            -1000024: "~chi1-", -1000023: "~chi20",
            -1000037: "~chi2-",- 1000025: "~chi30",
        }

        self.ids.update ( self.susyIDs )

    def rootColor( self, name ):
        """ find the default colors for <name>, ROOT version """
        from ROOT import kGreen,kOrange,kRed,kBlue,kBlack
        colors = { "orange": kOrange+3, "blue": kBlue+3, "red": kRed+2,
                   "black": kBlack, "green": kGreen+3 }
        c = self.namedColor ( name )
        return colors[c]

    def rgbColor( self, name, bold=False ):
        """ find the default colors for <name>, rgb version """
        boldcolors = { "orange": "#c56f18", "blue": "#000099", "red": "#990000",
                   "black": "#000000", "green": "#005500" }
        colors = { "orange": "#e58f38", "blue": "#0000ff", "red": "#ff0000",
                   "black": "#000000", "green": "#009900" }
        c = self.namedColor ( name )
        if bold:
            return boldcolors[c]
        return colors[c]

    def texColor( self, name ):
        """ find the default colors for <name>, rgb version """
        colors = { "orange": "{.4,.4,.1}", "blue": "{0,0,0.5}", "red": "{.5,0,0}",
                   "black": "{0,0,0}", "green": "{0,.5,0}" }
        c = self.namedColor ( name )
        # \color[rgb]{0,0,.5}\
        return f"\\color[rgb]{colors[c]}"

    def namedColor ( self, pid ):
        """ find the default colors for <name>, latex version 
        name can be an integer/pid, or a string/name. In case of a string,
        we will find the according pid.
        """
        if type(pid)==str:
            if not pid in self.names:
                print ( "[sparticleNames.texColor] %s not in names" % name )
                sys.exit(-1)
            pid = self.names[pid]
        pid = abs(pid)
        if pid in [ 1000001, 1000002, 1000003, 1000004, \
                    2000001, 2000002, 2000003, 2000004 ]:
            return "blue"
        if pid in [ 1000005, 1000006, 2000005, 2000006 ]:
            return "blue"
        if pid in [ 1000022, 1000023, 1000025, 1000035, \
                    1000024, 1000037 ]:
            return "green"
        if pid in [ 1000011, 1000012, 1000013, 1000014, \
                    1000015, 1000016, 2000011, 2000012, \
                    2000013, 2000014, 2000015, 2000016 ]:
            return "orange"
        return "black"

    def initXIDs ( self ):
        """ define the X notation """
        self.xIDs = { 1000001: "X_{d}", 2000001: "X_{d}^{2}",
            1000002: "X_{u}", 2000002: "X_{u}^{2}",
            1000003: "X_{s}", 2000003: "X_{s}^{2}",
            1000004: "X_{c}", 2000004: "X_{c}^{2}",
            1000005: "X_{b}", 2000005: "X_{b}^{2}",
            1000006: "X_{t}", 2000006: "X_{t}^{2}",
            1000011: "X_{e}", 1000012: "X_{\\nu}^{e}",
            1000013: "X_{\\mu}", 1000014: "X_{\\nu}^{\\mu}",
            1000015: "X_{\\tau}", 1000016: "X_{\\nu}^{\\tau}",
            2000011: "X_{e}", 2000012: "X_{\\nu}^{e}",
            2000013: "X_{\\mu}", 2000014: "X_{\\nu}^{\\mu}",
            2000015: "X_{\\tau}", 2000016: "X_{\\nu}^{\\tau}",
            1000021: "X_{g}", 1000022: "X_{Z}^{1}",
            1000024: "X_{W}", 1000023: "X_{Z}^{2}",
            1000037: "X_{W}^{2}", 1000025: "X_{Z}^{3}",
            1000035: "X_{Z}^{3}",
            -1000001: "#bar{X}_{d}", -2000001: "#bar{X}_{d}^{2}",
            -1000002: "#bar{X}_{u}", -2000002: "#bar{X}_{u}^{2}",
            -1000003: "#bar{X}_{s}", -2000003: "#bar{X}_{s}^{2}",
            -1000004: "#bar{X}_{c}", -2000004: "#bar{X}_{c}^{2}",
            -1000005: "#bar{X}_{b}", -2000005: "#bar{X}_{b}^{2}",
            -1000006: "#bar{X}_{t}", -2000006: "#bar{X}_{t}^{2}",
            -1000011: "#bar{X}_{e}", -1000012: "#bar{X}_{\\nu}^{e}",
            -1000013: "#bar{X}_{\\mu}", -1000014: "#bar{X}_{\\nu^{\\mu}",
            -1000015: "#bar{X}_{\\tau}", -1000016: "#bar{X}_{\\nu}^{\\tau}",
            -2000011: "#bar{X}_{e}", -2000012: "#bar{X}_{\\nu}^{e}",
            -2000013: "#bar{X}_{\\mu}", -2000014: "#bar{X}_{\\nu}^{\\mu}",
            -2000015: "#bar{X}_{\\tau}", -2000016: "#bar{X}_{\\nu}^{\\tau}",
            -1000021: "#bar{X}_{g}", -1000022: "#bar{X}_{Z}^{1}",
            -1000024: "#bar{X}_{W}", -1000023: "#bar{X}_{Z}^{2}",
            -1000037: "#bar{X}_{W}^{2}", -1000025: "~#bar{X}_{Z}^{3}",
            -1000035: "#bar{X}_{Z}^{3}"
        }
        self.ids.update ( self.xIDs )

    def __init__ ( self, susy=False ):
        """ Defines the ids and the names
        :param susy: use SUSY notation
        """
        self.susy = susy
        self.ids={
            1: "d", 2: "u", 3: "s", 4: "c", 5: "b", 6: "t", 11: "e^{-}", 13: "\\mu^{-}",
            15: "\\tau^{-}", 12: "\\nu", 14: "nu", 16:"nu", 21: "g", 22: "\\gamma",
            24: "W", 23:"Z", 25:"h1", 35: "h2", 36: "a0", 37: "h^{+}",
            -15: "\\tau^{+}", -13: "\\mu^{+}", -11: "e", -37: "h^{-}", -24: "W^{-}",
            -23: "Z", -25: "h1", -35: "h2", -36: "a0", -22: "\\gamma",
            -21: "g", -16: "\\nu", -14: "\\nu", -12: "nu", -1: "d",
            -2: "u", -3: "s", -4: "c", -5: "b", -6: "t" }
        if self.susy:
            self.initSUSY()
        else:
            self.initXIDs()

        ## make sure the inversions are defined as well
        self.names={}
        for (key,value) in self.ids.items():
            self.names[value]=key

    def isSM ( self, pid ):
        """ is pid a standard model pid? """
        if abs(pid)<100000:
            return True
        return False

    def rootify ( self, name ):
        """ rootify <name>, currently not doing anything """
        return name

    def htmlify ( self, name, addBrackets ):
        """ htmlify <name> """
        import re
        html = name

        replacements = { "\\mu": "&mu;", "\\nu": "&nu;", "\\tau": "&tau;", 
                         "\\gamma": "&gamma;" }
        for k,v in replacements.items():
            html = html.replace(k,v)

        while True:
            m = re.search ( "_{[A-Za-z&;]*}", html)
            if m == None:
                break
            repl = html[m.start()+2:m.end()-1]
            html = html[:m.start()]+"<sub>"+repl+"</sub>"+html[m.end():]

        while True:
            m = re.search ( "\^{[0-9&;A-Za-z-+]*}", html)
            if m == None:
                break
            repl = html[m.start()+2:m.end()-1]
            html = html[:m.start()]+"<sup>"+repl+"</sup>"+html[m.end():]

        if addBrackets:
            html = f"({html})"
        html = html.replace("#bar{X}","x&#772;" )
        html = html.replace("~","")
        html = html.replace("_{&nu;","<sub>&nu;</sub>")
        if html == "nu":
            html="&nu;"
        # print ( "htmlify", name, "->", html )

        return html

    def rootName ( self, pid, addSign = False ):
        """ format the name for ROOT """
        return self.rootify ( self.name ( pid, addSign ) )

    def htmlName ( self, pid, addSign=False, addBrackets = False ):
        """ format the name for html """
        return self.htmlify ( self.name ( pid, addSign ), addBrackets )

    def texName ( self, pid, addSign=False, addDollars = False, addBrackets = False ):
        """ format the name for tex """
        n = self.name ( pid, addSign )
        n = n.replace ( "#", "\\" )
        if addSign:
            if "+-" in n:
                n = n.replace("+-","" )
                n += "^{\\pm}"
        else:
            n = n.replace("^{-}","").replace("^{+}","")
        if addDollars:
            n = "$" + n + "$"
        if addBrackets:
            n = "(" + n + ")"
        return n

    def name ( self, pid, addSign=False ):
        """ get the name for a particle id """
        if type(pid) in [ tuple, set, list ]:
            ret=[]
            for p in pid:
                ret.append ( self.name ( p, addSign ) )
            return ", ".join ( ret )

        if type(pid) == str and pid.startswith("+-"):
            n = self.name(int(pid[2:]))
            return "+-"+n
        if not pid in self.ids and not abs(pid) in self.ids:
            return str(pid)
        if not pid in self.ids:
            return self.ids[abs(pid)]
        return self.ids[pid]

    def asciiName ( self, pid ):
        """ get the ascii version of the name """
        ret = self.name ( pid )
        ret = ret.replace("_","").replace("{","").replace("}","").replace("^","").replace("\\","")
        return ret

    def pid ( self, name ):
        """ get the pid for a particle name """
        if not name in self.names:
            return 0
        return self.names[name]

    def has ( self, i ):
        """ do we have particle? can be pid or name """
        if i in self.names: return True
        if i in self.ids: return True
        return False

    def particleType ( self, pid ):
        """ categorizes sparticles """
        q=abs(pid)
        if q>1000000 and q<1000005:
            return "q"
        if q>2000000 and q<2000005:
            return "q"
        if q in [ 1000005, 2000005 ]:
            return "b"
        if q in [ 1000006, 2000006 ]:
            return "t"
        if q==1000021:
            return "g"
        if q in [ 1000022, 1000023, 1000025, 1000035, 1000024, 1000037 ]:
            return "n"
        if q in [ 1000011, 1000013, 1000015, 2000011, 2000013, 2000015 ]:
            return "l"
        if q in [ 1000012, 1000014, 1000016, 2000012, 2000014, 2000016 ]:
            return "l"
        return str(q)

    def shortName ( self, productiontuple ):
        """ assign a particle category to a tuple of two particle pids """
        p1,p2=abs( productiontuple[0] ),abs( productiontuple[1] )
        # p1,p2= productiontuple
        q1,q2=self.particleType ( p1 ), self.particleType ( p2 )
        if q1>q2: q1,q2=q2,q1 ## swap, give a canonical order
        return q1+q2

    def longName ( self, letter ):
        """ gives long names to particle categories """
        if letter=="l": return "slepton"
        if letter=="n": return "weakino"
        if letter=="q": return "squark"
        if letter=="t": return "stop"
        if letter=="b": return "sbottom"
        if letter=="g": return "gluino"
        return "?"

    def tilde ( self, text ):
        """ put a tilde over text """
        return "<math display='inline'><mover><mi>%s</mi><mo stretchy='true'>~</mo></mover></math>" % text

    def sub ( self, text ):
        return "<sub>"+text+"</sub>"

    def sup ( self, text ):
        return "<sup>"+text+"</sup>"

    """
    def toHtml ( self, name ):
        # translate particle names to html code
        #if name=="~chi2+":
        #    name=self.tilde("&chi;")+"xxx" ## sup("+")+"2"#+sub("2")
        #if name=="~chi1+": return self.tilde("&chi;")+"1+" ## sup("+")+"2"#+sub("2")
        #if name=="~chi30": return self.tilde("&chi;")+sup("0")#+sub("2")
        name=name.replace("_eL","<sub>eL</sub>")
        name=name.replace("_muL","<sub>muL</sub>")
        name=name.replace("_tauL","<sub>tauL</sub>")
        name=name.replace("chi10","chi<sub>1</sub><sup>0</sup>")
        name=name.replace("chi20","chi<sub>2</sub><sup>0</sup>")
        name=name.replace("chi30","chi<sub>3</sub><sup>0</sup>")
        name=name.replace("chi40","chi<sub>4</sub><sup>0</sup>")
        name=name.replace("chi50","chi<sub>5</sub><sup>0</sup>")
        name=name.replace("chi1+","chi<sub>1</sub><sup>+</sup>")
        name=name.replace("chi2+","chi<sub>2</sub><sup>+</sup>")
        name=name.replace("chi3+","chi<sub>3</sub><sup>+</sup>")
        name=name.replace("chi1-","chi<sub>1</sub><sup>-</sup>")
        name=name.replace("chi2-","chi<sub>2</sub><sup>-</sup>")
        name=name.replace("chi3-","chi<sub>3</sub><sup>-</sup>")
        name=name.replace("chi","&chi;")
        name=name.replace("nu","&nu;")
        name=name.replace("mu","&mu;")
        name=name.replace("tau","&tau;")
        name=name.replace("_L","<sub>L</sub>")
        name=name.replace("uL","u<sub>L</sub>")
        name=name.replace("dL","d<sub>L</sub>")
        name=name.replace("cL","c<sub>L</sub>")
        name=name.replace("sL","s<sub>L</sub>")
        name=name.replace("uR","u<sub>R</sub>")
        name=name.replace("dR","d<sub>R</sub>")
        name=name.replace("cR","c<sub>R</sub>")
        name=name.replace("sR","s<sub>R</sub>")
        name=name.replace("_R","<sub>R</sub>")
        name=name.replace("_1","<sub>1</sub>")
        name=name.replace("A0","A<sup>0</sup>")
        name=name.replace("H+","H<sup>+</sup>")
        name=name.replace("_2","<sub>2</sub>")
        name=name.replace("b1","b<sub>1</sub>")
        name=name.replace("b2","b<sub>2</sub>")
        name=name.replace("t1","t<sub>1</sub>")
        name=name.replace("t2","t<sub>2</sub>")
        name=name.replace("^*","<sup>*</sup>")
        #name=name.replace("+","<sup>+</sup>")
        if name.find("~")==0:
            if name.find("<su")==-1:
                name=self.tilde(name[1:])
            else:
                pos=name.find("<su")
                name=self.tilde(name[1:pos])+name[pos:]
        # print name,"<br>"
        return "<nobr>"+name+"</nobr>"
    """


if __name__ == "__main__":
    """ as a script, we simply print out the paths """
    print ( "sparticle names" )
    namer = SParticleNames()
    ctr=0
    f=open("index.html","wt" )
    f.write ( "<html><body>\n" )
    for (key,value) in namer.ids.items():
       ctr+=1
       print ()
       line = "pid %d: %s" % ( key, namer.htmlName ( key ) )
       print ( line )
       print ()
       f.write ( line + "<br>\n" )
    f.close()
       # print ( "%8d %8s   |" % (key,value), end="" )
       #if ctr==3:
       #  print ()
       #  ctr=0
