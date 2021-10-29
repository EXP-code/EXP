"""
psp_io.py
   driver for input classes and handling particle data in Python

MSP 25 Oct 2021 Overhauled version


TODO

"""
import numpy as np

try:
    # requires yaml support: likely needs to be installed.
    import yaml
except ImportError:
    raise ImportError("You will need to 'pip install pyyaml' to use this reader.")
    


class Input():
    """Input class to adaptively handle various EXP output formats.

    inputs
    ---------------
    filename : str
        the input filename to be read
    comp     : str, optional
        the name of the component for which to extract data. If None, will read primary header and exit.
    legacy   : bool, default=True
        if True, return attributes rather than a dictionary of particle data
    verbose  : int, default 0
        verbosity flag.
    
    returns
    ---------------
    self        : Input instance
      .header   : dict, all header values pulled from the file
        the .keys() are the names of each component
        each component has a dictionary of values, including 'parameters'
        the details of the force calculation are in 'force'
      .filename : str, the filename that was read
      .comp     : str, name of the component
      .time     : float, the time in the output file
      .data     : dictionary, with keys:
        x       : float, the x position
        y       : float, the y position
        z       : float, the z position
        vx      : float, the x velocity
        vy      : float, the y velocity
        vz      : float, the z velocity
        mass    : float, the mass of the particle
        index   : int, the integer index of the particle
        potE    : float, the potential energy value

      ---or, if legacy=True---
      
      .xpos     : float, the x position
      .ypos     : float, the y position
      .zpos     : float, the z position
      .xvel     : float, the x velocity
      .yvel     : float, the y velocity
      .zvel     : float, the z velocity
      .mass     : float, the mass of the particle
      .indx     : int, the integer index of the particle
      .pote     : float, the potential energy value


    """
    def __init__(self, filename,comp=None, legacy=True,verbose=0):
        """the main driver for the class"""

        # auto-determine the file type
        if 'SPL.' in filename:
            self.style = 'SPL'
        elif 'OUT.' in filename:
            self.style = 'OUT'
        else:
            self.style = 'unknown'
        
        if self.style=='SPL':
            I = SPLInput(filename,comp=comp,verbose=verbose)
        elif self.style=='OUT':
            I = OUTInput(filename,comp=comp,verbose=verbose)
        else:
            raise ValueError('File type not supported for file "{}"'.format(filename))

        # expose the header
        self.header = I.header
        
        
        # what is the ideal legacy error handling?
        if I.comp==None:
            return
            
        if legacy:
            self.mass = I.data['m']
            self.xpos = I.data['x']
            self.ypos = I.data['y']
            self.zpos = I.data['z']
            self.xvel = I.data['vx']
            self.yvel = I.data['vy']
            self.zvel = I.data['vz']
            self.pote = I.data['potE']
            
            if (I.header['dark']['parameters']['indexing']):
                self.indx = I.data['index']

        else:
            self.data = I.data

        self.filename = I.filename
        self.comp     = I.comp
        self.time     = I.time
        

        

class OUTInput:
    """Input class to adaptively handle OUT. format specifically

    inputs
    ---------------
    filename : str
        the input filename to be read
    comp     : str, optional
        the name of the component for which to extract data. If None, will read primary header and exit.
    verbose  : int, default 0
        verbosity flag.
    
    returns
    ---------------
    self        : Input instance
      .header   : dict, all header values pulled from the file
        the .keys() are the names of each component
        each component has a dictionary of values, including 'parameters'
        the details of the force calculation are in 'force'
      .filename : str, the filename that was read
      .comp     : str, name of the component
      .time     : float, the time in the output file
      .data     : dictionary, with keys:
        x       : float, the x position
        y       : float, the y position
        z       : float, the z position
        vx      : float, the x velocity
        vy      : float, the y velocity
        vz      : float, the z velocity
        mass    : float, the mass of the particle
        index   : int, the integer index of the particle
        potE    : float, the potential energy value

    """
    def __init__(self, filename,comp=None,verbose=0):
        """the main driver"""
        self.verbose  = verbose       
        self.filename = filename
        
        # initial check for file validity
        try:
            self.f = open(self.filename, 'rb')
        except Exception:
            raise IOError('Failed to open "{}"'.format(filename))

        # do an initial read of the header
        self.primary_header = dict()

        # initialise dictionaries
        self.comp_map       = dict()
        self.header      = dict()
        
        self._read_primary_header()

        self.comp = comp
        _comps = list(self.header.keys())
 
        # if a component is defined, retrieve data
        if comp != None:
            if comp not in _comps:
                raise IOError('The specified component does not exist.')
                
            else:
                self.comp = comp
                self.data = self._read_component_data(self.filename,
                                                      self.header[self.comp]['nbodies'],
                                                      int(self.header[self.comp]['data_start']))
        # wrapup
        self.f.close()
        
    def _read_primary_header(self):
        """read the primary header from an OUT. file"""

        self._check_magic_number()

        # reset to beginning and read current time
        self.f.seek(0)
        self.time, = np.fromfile(self.f, dtype='<f8', count=1)
        self._nbodies_tot, self._ncomp = np.fromfile(self.f, dtype=np.uint32,count=2)

        data_start = 16
        
        for comp in range(0,self._ncomp):
            self.f.seek(data_start)
            next_comp = self._read_out_component_header()
            data_start = next_comp

    def _summarise_primary_header(self):
        """a short summary of what is in the file"""

        ncomponents = len(self.header.keys())
        comp_list   = list(self.header.keys())
        print("Found {} components.".format(ncomponents))

        for n in range(0,ncomponents):
            print("Component {}: {}".format(n,comp_list[n]))
            

    def _read_out_component_header(self):
        """read in the header for a single component, from an OUT. file"""
        
        #_ = f.tell()  # byte position of this component

        if self._float_len == 4:
            _1,_2, nbodies, nint_attr, nfloat_attr, infostringlen = np.fromfile(self.f, dtype=np.uint32, count=6)            
        else:
            nbodies, nint_attr, nfloat_attr, infostringlen = np.fromfile(self.f, dtype=np.uint32, count=4)

        # information string from the header
        head = np.fromfile(self.f, dtype=np.dtype((np.bytes_, infostringlen)),count=1)
        head_normal = head[0].decode()
        head_dict = yaml.safe_load(head_normal)
        
        # deprecated backward compatibility here: see frozen versions if this is an old file

        comp_data_pos = self.f.tell()  # byte position where component data begins

        # the default fields are (m, x, y, z, vx, vy, vz, p)
        nfields = 8
        comp_length = nbodies * (8 * int(head_dict['parameters']['indexing']) +
                                 self._float_len * nfields +
                                 4 * nint_attr +
                                 self._float_len * nfloat_attr)
        comp_data_end = self.f.tell() + comp_length  # byte pos. of comp. data end

        head_dict['nint_attr']   = nint_attr
        head_dict['nfloat_attr'] = nfloat_attr
        head_dict['nbodies']     = nbodies
        head_dict['data_start']  = comp_data_pos
        head_dict['data_end']    = comp_data_end

        self.header[head_dict['name']] = head_dict

        # specifically look for indexing
        try:
            self.indexing = head_dict['parameters']['indexing']
        except:
            self.indexing = head_dict['indexing']=='true'

        return comp_data_end
                
    
    def _check_magic_number(self):
        """check the magic number to see if a file is float or double"""

        self.f.seek(16)  # find magic number
        cmagic, = np.fromfile(self.f, dtype=np.uint32, count=1)

        # check if it is float vs. double
        if cmagic == 2915019716:
            self._float_len = 4
            self._float_str = 'f'
        else:
            self._float_len = 8
            self._float_str = 'd'
        


    def _read_component_data(self,filename,nbodies,offset):
        """read in all data for component"""

        dtype_str = []
        colnames = []
        if self.header[self.comp]['parameters']['indexing']:
            # if indexing is on, the 0th column is Long
            dtype_str = dtype_str + ['l']
            colnames  = colnames + ['index']

        dtype_str = dtype_str + [self._float_str] * 8
        colnames = colnames + ['m', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'potE']

        dtype_str = dtype_str + ['i'] * self.header[self.comp]['nint_attr']
        colnames = colnames + ['i_attr{}'.format(i)
                               for i in range(self.header[self.comp]['nint_attr'])]

        dtype_str = dtype_str + [self._float_str] * self.header[self.comp]['nfloat_attr']
        colnames = colnames + ['f_attr{}'.format(i)
                               for i in range(self.header[self.comp]['nfloat_attr'])]

        dtype = np.dtype(','.join(dtype_str))

        out = np.memmap(filename,
                        dtype=dtype,
                        shape=(1, nbodies),
                        offset=offset,
                        order='F', mode='r')

        tbl = dict()
        for i, name in enumerate(colnames):
            tbl[name] = np.array(out['f{}'.format(i)][0], copy=True)

        del out  # close the memmap instance

        return tbl




class SPLInput:
    """Input class to adaptively handle SPL. format specifically

    inputs
    ---------------
    filename : str
        the input filename to be read
    comp     : str, optional
        the name of the component for which to extract data. If None, will read primary header and exit.
    verbose  : int, default 0
        verbosity flag.
    
    returns
    ---------------
    self        : Input instance
      .header   : dict, all header values pulled from the file
        the .keys() are the names of each component
        each component has a dictionary of values, including 'parameters'
        the details of the force calculation are in 'force'
      .filename : str, the filename that was read
      .comp     : str, name of the component
      .time     : float, the time in the output file
      .data     : dictionary, with keys:
        x       : float, the x position
        y       : float, the y position
        z       : float, the z position
        vx      : float, the x velocity
        vy      : float, the y velocity
        vz      : float, the z velocity
        mass    : float, the mass of the particle
        index   : int, the integer index of the particle
        potE    : float, the potential energy value

    """
    def __init__(self, filename,comp=None,verbose=0):
        """the main driver, see above for parameters"""
        self.verbose  = verbose
        self.filename = filename

        # initial check for file validity
        try:
            self.f = open(self.filename, 'rb')
        except Exception:
            raise IOError('Failed to open "{}"'.format(filename))

        # do an initial read of the header
        self.primary_header = dict()

        # initialise dictionaries
        self.component_map = dict()
        self.header        = dict()
        
        self._read_primary_header()

        self.comp = comp
        _comps = list(self.header.keys())

        if comp == None:
            self._summarise_primary_header()
            return
        
        # if a component is defined, retrieve data
        if comp != None:
            if comp not in _comps:
                raise IOError('The specified component does not exist.')
                
            else:
                self.comp = comp
                
                # set up assuming the directory of the main file is the same
                # as the subfiles. could add verbose flag to warn?
                self.indir = filename.split('SPL')[0]

                # now we can query out a specific component
                self._make_spl_file_list(self.comp)

                # given the comp, pull the data.
                self._read_spl_component_data()

        # wrapup
        self.f.close()


    def _read_primary_header(self):
        """read the primary header from an SPL. file"""

        self._check_magic_number()
        
        self.f.seek(0)
        self.time, = np.fromfile(self.f, dtype='<f8', count=1)
        self._nbodies_tot, self._ncomp = np.fromfile(self.f, dtype=np.uint32,count=2)

        data_start = 16 # guaranteed first component location...

        # now read the component headers
        for comp in range(0,self._ncomp):
            self.f.seek(data_start) 
            next_comp = self._read_spl_component_header()               
            data_start = next_comp 
    

    def _read_spl_component_header(self):
        """read in the header for a single component, from an SPL. file"""

        data_start = self.f.tell()
        # manually do headers
        _1,_2,self.nprocs, nbodies, nint_attr, nfloat_attr, infostringlen = np.fromfile(self.f, dtype=np.uint32, count=7)

        # need to figure out what to do with nprocs...it has to always be the same, right?
        head = np.fromfile(self.f, dtype=np.dtype((np.bytes_, infostringlen)),count=1)
        head_normal = head[0].decode()
        head_dict = yaml.safe_load(head_normal)

        
        head_dict['nint_attr']   = nint_attr
        head_dict['nfloat_attr'] = nfloat_attr
        # data starts at ...
        next_comp = 4*7 + infostringlen + data_start
            
        self.component_map[head_dict['name']]  = next_comp
        self.header[head_dict['name']]         = head_dict

        next_comp +=  self.nprocs*1024

        # specifically look for indexing
        try:
            self.indexing = head_dict['parameters']['indexing']
        except:
            self.indexing = head_dict['indexing']=='true'
        
        return next_comp

                
    def _make_spl_file_list(self,comp):
        
        self.f.seek(self.component_map[comp])
        
        PBUF_SZ = 1024
        PBUF_SM = 32

        self.subfiles = []
    
        for procnum in range(0,self.nprocs):
            PBUF = np.fromfile(self.f, dtype=np.dtype((np.bytes_, PBUF_SZ)),count=1)
            subfile = PBUF[0].split(b'\x00')[0].decode()
            self.subfiles.append(subfile)

    def _summarise_primary_header(self):
        """a short summary of what is in the file"""

        ncomponents = len(self.header.keys())
        comp_list   = list(self.header.keys())
        print("Found {} components.".format(ncomponents))

        for n in range(0,ncomponents):
            print("Component {}: {}".format(n,comp_list[n]))
            

    def _read_spl_component_data(self):
        
        FullParticles = dict()

        # first pass: get everything into memory
        for n in range(0,len(self.subfiles)):
            FullParticles[n] = dict()

            if self.verbose>1:
                print('spl_io._read_spl_component_data: On file {} of {}.'.format(n,len(self.subfiles)))
            
            tbl = self._handle_spl_subfile(self.subfiles[n])

            for k in tbl.keys():
                FullParticles[n][k] = tbl[k]

        # construct a single dictionary for the particles
        self.data = dict()
        for k in tbl.keys():
            self.data[k] = np.concatenate([FullParticles[n][k] for n in range(0,len(self.subfiles))])

        # cleanup...
        del FullParticles
        
            
    def _check_magic_number(self):

        self.f.seek(16)  # find magic number
        cmagic, = np.fromfile(self.f, dtype=np.uint32, count=1)

        # check if it is float vs. double
        if cmagic == 2915019716:
            self._float_len = 4
            self._float_str = 'f'
        else:
            self._float_len = 8
            self._float_str = 'd'
        

    def _handle_spl_subfile(self,subfilename):

        subfile = open(self.indir+subfilename,'rb')
        nbodies, = np.fromfile(subfile, dtype=np.uint32, count=1)
        subfile.close() # close the opened file

        # read the data
        tbl = self._read_component_data(self.indir+subfilename,nbodies,4)
        # the offset is always fixed in SPL subfiles

        return tbl


    def _read_component_data(self,filename,nbodies,offset):
        """read in all data for component"""
        
        dtype_str = []
        colnames  = []
        if self.header[self.comp]['parameters']['indexing']:
            # if indexing is on, the 0th column is Long
            dtype_str = dtype_str + ['l']
            colnames  = colnames + ['index']
        
        dtype_str = dtype_str + [self._float_str] * 8
        colnames = colnames + ['m', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'potE']
        
        dtype_str = dtype_str + ['i'] * self.header[self.comp]['nint_attr']
        colnames = colnames + ['i_attr{}'.format(i)
                               for i in range(self.header[self.comp]['nint_attr'])]
        
        dtype_str = dtype_str + [self._float_str] * self.header[self.comp]['nfloat_attr']
        colnames = colnames + ['f_attr{}'.format(i)
                               for i in range(self.header[self.comp]['nfloat_attr'])]
        
        dtype = np.dtype(','.join(dtype_str))
        
        out = np.memmap(filename,
                        dtype=dtype,
                        shape=(1, nbodies),
                        offset=offset,
                        order='F', mode='r')
        
        tbl = dict()
        for i, name in enumerate(colnames):
            tbl[name] = np.array(out['f{}'.format(i)][0], copy=True)
        
        del out  # close the memmap instance
        
        return tbl


