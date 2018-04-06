#!/bin/sh
platex kh_inst.tex
dvipdfmx kh_inst.dvi
open kh_inst.pdf