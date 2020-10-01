This folder is all about combining different analyses, and finding the best
combination. Including the artificial model builder.

Random walker
^^^^^^^^^^^^^

  * Randomly builds models and compares them against the SMS results
  * Try: ./walker.py -h

Accelerator
^^^^^^^^^^^

  * The code that trains neural networks to predict the Z, plus code for the gradient ascent.

  * Called from the command line, can be used to train the Z score prediction network
  * Try: ./accelerator.py -h

Hiscore keeper
^^^^^^^^^^^^^^
  * The hiscore.py script can be used to collect the hiscores from the different walkers, 
    trim them, compute analysis contributions.  
  * Try: ./hiscore.py -h
  * Use plotHiscore.py to produce web pages like http://www.hephy.at/user/wwaltenberger/models/

Scanner
^^^^^^^

  * Performs parameter scans around the winning protomodel

Llhdscanner
^^^^^^^^^^^

  * Code that performs likelihood scans per analysis and topology, to see if
    it really looks like a signal
  * Llhdplot plots the scans

History
^^^^^^^

 * Used for plotting a specific random walk

updateHiscores
^^^^^^^^^^^^^^

etc etc
