#!/usr/bin/env python3

# mostly taken from http://code.google.com/p/latexmath2png/

import os
import re
import sys
import shutil
import tempfile
import subprocess

class Latex(object):
  INLINE = r'''
\documentclass[fleqn]{article}
\usepackage{amssymb,amsmath,bm,color}
\usepackage[latin1]{inputenc}
\usepackage[active,textmath]{preview}
\begin{document}
\thispagestyle{empty}
\mathindent0cm
\parindent0cm
%s
\end{document}
'''
  BLOCK = r'''
\documentclass[fleqn]{article}
\usepackage{amssymb,amsmath,bm,color}
\usepackage[latin1]{inputenc}
\begin{document}
\thispagestyle{empty}
\mathindent0cm
\parindent0cm
%s
\end{document}
'''

  def __init__(self, doc, dpi=120):
    self.doc = doc
    self.dpi = dpi

  def write(self):
    inline = bool(re.match(r'^\$[^$]*\$$', self.doc))
    if inline:
      TEX = self.INLINE
    else:
      TEX = self.BLOCK
    if not shutil.which('dvipng'):
      TEX = self._prepare_pdf_template(TEX, inline)

    try:
      workdir = tempfile.gettempdir()
      fd, texfile = tempfile.mkstemp('.tex', 'eq', workdir, True)

      with os.fdopen(fd, 'w+') as f:
        f.write(TEX % self.doc)

      png, depth = self.convert_file(texfile, workdir)
      if not inline:
        depth = 0
      return png, depth

    finally:
      if os.path.exists(texfile):
        os.remove(texfile)

  def _prepare_pdf_template(self, TEX, inline):
    if inline:
      return TEX.replace('[active,textmath]', '[active,tightpage,textmath]')
    TEX = TEX.replace(r'\usepackage[latin1]{inputenc}',
        r'\usepackage[latin1]{inputenc}' + '\n' + r'\usepackage[active,tightpage]{preview}')
    return TEX.replace('%s', r'\begin{preview}%s\end{preview}')

  def convert_file(self, infile, workdir):

    try:
      dvipng = shutil.which('dvipng')
      if dvipng:
        cmd = f'latex -halt-on-error -output-directory {workdir} {infile}'
      else:
        cmd = f'pdflatex -halt-on-error -interaction=nonstopmode -output-directory {workdir} {infile}'

      p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      sout, serr = p.communicate()
      # Something bad happened, abort
      if p.returncode != 0:
          raise Exception('latex error' if dvipng else 'pdflatex error', serr, sout)

      # Convert the DVI or PDF file to PNG's
      if dvipng:
        dvifile = infile.replace('.tex', '.dvi')
        pngfile, depth = self._convert_with_dvipng(dvifile, workdir)
      else:
        pdffile = infile.replace('.tex', '.pdf')
        pngfile, depth = self._convert_with_poppler(pdffile)

      png = open(pngfile,'rb').read()
      return png, depth

    finally:
      # Cleanup temporaries
      basefile = infile.replace('.tex', '')
      tempext = '.aux', '.dvi', '.log', '.pdf', '.png'
      for te in tempext:
        tempfile = basefile + te
        if os.path.exists(tempfile):
          os.remove(tempfile)

  def _convert_with_dvipng(self, dvifile, workdir):
    pngfile = os.path.join(workdir, os.path.basename(dvifile).replace('.dvi', '.png'))
    cmd = f"dvipng -T tight --depth -D {int(self.dpi)} -z 9 -bg White -o {pngfile} {dvifile}"
    #cmd = "dvipng -T tight --depth -D %i -z 9 -bg Transparent -o %s %s" % (self.dpi, pngfile, dvifile)
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    sout, serr = p.communicate()
    if p.returncode != 0:
      raise Exception('dvipng error', serr)

    try:
      depth = int(re.search(r'\[1 depth=(-?\d+)\]', sout).group(1))
    except:
      depth = 0

    return pngfile, depth

  def _convert_with_poppler(self, pdffile):
    if not shutil.which('pdftoppm'):
      raise Exception('dvipng not found and pdftoppm is unavailable; install poppler')

    pngbase = pdffile.replace('.pdf', '')

    cmd = f"pdftoppm -png -r {int(self.dpi)} -singlefile {pdffile} {pngbase}"
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    sout, serr = p.communicate()
    if p.returncode != 0:
      raise Exception('pdftoppm error', serr)

    return pngbase + '.png', 0


__cache = {}

def tex2png(eq, **kwargs):
  if not eq in __cache:
    __cache[eq] = Latex(eq, **kwargs).write()
  return __cache[eq]

if __name__ == '__main__':
  src = sys.argv[1]
  print ( f'Equation is: {src}' )
  p = Latex(src,600).write()
  f= open ( "out.png", "wb" )
  f.write ( p[0] )
  f.close()
