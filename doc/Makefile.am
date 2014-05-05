## Makefile for EXP docs
##
## Author: Martin Weinberg on 03/03/02
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

## Process this file with automake to produce Makefile.in

DOCS = *.doc *.html

FORMAT_STR = "s/\(PROJECT_NUMBER[ ]*\=\).*/\1 $(PACKAGE_VERSION)/g"

all: docs-html

if DOXYGEN_FOUND
docs-html:	$(DOCS)	
		cp exp.cfg tmp.cfg
		sed -e $(FORMAT_STR) < tmp.cfg > exp.cfg
		$(DOXYGEN) exp.cfg
endif

docs-clean:
	rm -Rf latex html


