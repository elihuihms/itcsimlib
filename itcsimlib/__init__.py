__version__		= "0.6.0"
__copyright__	= """Copyright (C) Elihu Ihms

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA."""

__license__	= u'GPL v2'
__author__	= u'Elihu Ihms'
__author_email__	= u'mail@elihuihms.com'

MATPLOTLIB_BACKEND = None

from .itc_sim		import ITCSim
from .itc_model	import ITCModel
from .itc_fit		import ITCFit
from .itc_calc		import ITCCalc
from .itc_grid		import ITCGrid
