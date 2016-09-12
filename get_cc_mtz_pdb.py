import iotbx
from iotbx import pdb
import libtbx.phil
from cctbx.array_family import flex
import math
import libtbx.phil.command_line
from libtbx import runtime_utils
from libtbx.str_utils import format_value
from libtbx.utils import Sorry
from phenix.autosol.delete_file import delete_file
from phenix.autosol.copy_to_temp_dir import copy_to_temp_dir
from phenix.autosol.copy_from_temp_dir import copy_from_temp_dir
from phenix.autosol.trim_file_name import trim_file_name
from phenix.autosol.run_resolve import run_resolve
import sys, os, re, string,time
from phenix.autosol.UserMethods import GeneralMethods
from phenix.utilities.file_name_methods import file_name_methods
from phenix.utilities.composite_params import get_composite_master_params
from phenix.utilities.is_debug import is_debug
from phenix.utilities.catenate_equals import catenate_equals
from mmtbx import pdbtools

master_params="""
get_cc_mtz_pdb {

      pdb_in = None
        .type = path
        .help = '''PDB file with coordinates to evaluate'''
        .short_caption = PDB file
        .style = bold noauto file_type:pdb

      atom_selection = None
      .type = str
      .help = "Any selection specified with atom_selection is applied to "
              "input model (pdb_in) before using the model."
              "NOTE: this option can result in confusing output because "
              "the model density near an atom and therefore the local CC value"
              "depends not just on the atoms being considered but"
              "also nearby atoms.  Therefore if you remove some atoms next "
              "to a residue with atom_selection then the CC for that residue"
              "may change."
      .input_size = 400
      .style = bold

      mtz_in = None
        .type = path
        .help = '''MTZ file with coefficients for a map '''
        .short_caption = Map coefficients
        .style = bold noauto file_type:mtz process_hkl child:map_labels:labin \
          child:d_min:resolution

      labin = ""
        .type = str
        .help = ''' Labin line for MTZ file with map coefficients.
                 This is optional if get_cc_mtz_pdb
                 can guess the correct coefficients
                 for FP PHI and FOM.  Otherwise specify:
                 LABIN FP=myFP PHIB=myPHI FOM=myFOM
                 where myFP is your column label for FP '''
        .short_caption = Data labels
        .style = bold noauto renderer:draw_map_coeffs_widget \
          parent:file_name:mtz_in OnChange:auto_update_label_choice \
          child:d_min:resolution

      offset_pdb = "offset.pdb"
        .type = path
        .help = "Output version of pdb file, offset to maximize "
                "correlation with mtz file"
        .style = bold

      resolution = 0.
        .type = float
        .help = "high-resolution limit for map calculation"
        .short_caption = High resolution
        .style = bold resolution

      use_only_refl_present_in_mtz = False
        .type = bool
        .help = ''' You can specify that only reflections present in
          your mtz file are used in the comparison.
          '''
        .short_caption = Use only reflections present in MTZ file

      scale = False
        .type = bool
        .help = '''If you set scale=True then get_cc_mtz_pdb applies
           an overall B factor and a delta_b for each atom beyond CB.
           '''
        .short_caption = Scale sidechain B-factors

      select_by_b = False
        .type = bool
        .help = '''select atoms by their b factors'''
        .short_caption = Select atoms by B-factors

      split_conformers = False
        .type = bool
        .help = '''If you want to have A and B conformers analyzed separately
               you can say split_conformers=True. Note that this requires that
               the entire A conformer for a residue must be before (or after)
               the entire B conformer for that residue.'''
        .short_caption = Split conformers

      fix_xyz = False
        .type = bool
        .help = '''If you want your PDB file compared to the map from
               your mtz file with no offsets at all (fixed position)
               then specify fix_xyz=True'''
        .short_caption = Fix PDB position

      fix_rad_max = False
        .type = bool
        .help = "If you want to use a fixed radius around all atoms "
                "for calculation of the correlation use fix_rad=True."
                "To set the value, use rad_max=xxx, otherwise it is set"
                "automatically.  If rad_max is set, fix_rad_max is set"
                "to True"
        .short_caption = Fix CC radius rad_max
        .style = hidden

      rad_max = None
        .type = float
        .help = "If you want to use a fixed radius around all atoms "
                "for calculation of the correlation use fix_rad=True."
                "To set the value, use rad_max=2.5, otherwise it is set"
                "automatically"
        .short_caption = CC radius rad_max
        .style = hidden

      any_offset = False
        .type = bool
        .help = '''You can search for a match with any offset even
          though this is not allowed by space-group symmetry'''
        .short_caption = Search for match with any offset

      chain_type = *PROTEIN DNA RNA
        .type = choice
        .help = "Chain type (for identifying main-chain and side-chain atoms)"
        .short_caption = Chain type
        .style = bold

      include scope phenix.utilities.common_params.gui_directory_params_str

      verbose = True
        .type = bool
        .help = '''Verbose output'''
        .short_caption = Verbose output

      quick   = False
        .type = bool
        .help = '''Skip the residue-by=residue correlations for a quick run'''
        .short_caption = Skip residue-by-residue correlations

      raise_sorry = False
        .type = bool
        .help = "Raise sorry if problems"
        .short_caption = raise sorry

      debug = False
        .type = bool
        .help = "Debugging output"
        .short_caption = Debugging output

      dry_run = False
        .type = bool
        .help = '''Just read in and check parameter names'''

      include scope libtbx.phil.interface.tracking_params

}
"""

class get_cc_mtz_pdb(GeneralMethods):
  def __init__(self,args,quiet=False,out=sys.stdout):
    command_name="get_cc_mtz_pdb"
    self.Name='get_cc_mtz_pdb'
    self.found_overall=0
    self.found_region=0
    from phenix.utilities import citations
    citations.add_citation('phenix','get_cc_mtz_pdb')
#to process pdb files
    command_line_interpreter = pdbtools.interpreter(command_name  = command_name,
                                         args          = args,
                                         log           = out,
                                         get_cc        = True)
    command_line_interpreter.set_ppf()
    xray_structure = command_line_interpreter.processed_pdb_file.xray_structure(
    show_summary = False)
    if(xray_structure is None):
      raise Sorry("Cannot extract xray_structure.")
    xs_a = xray_structure
    u_a     = xs_a    .extract_u_iso_or_u_equiv()
    eps = math.pi**2*8

    def calc_std(values, eps = None):
      if(eps is not None): values = values * eps
      values = values - flex.mean_default(values, None)
      values = flex.pow(values, 2)
      return math.sqrt(flex.sum(values)) / values.size()

    def select_std(values, mean, sigma, eps = None):
      if(eps is not None): values = values * eps
      return values > mean + sigma

    std = calc_std(u_a, eps)
    high_u_a = select_std(u_a, flex.mean_default(u_a, None), 3*std, eps)

    args=catenate_equals(args).new_args()
    args=self.special_cases(args)

    master_params=get_composite_master_params(
         method_list=['get_cc_mtz_pdb'],
         location_list=['phenix.command_line'])

    args=self.get_keyword_table(args,out=out)       # set self.keyword_table

    summary,header=self.get_summary_and_header(command_name)
    done,master_params,params,changed_params,help=self.get_params(
        command_name,master_params,args,out=out)
    if done: return
    if not quiet: print >>out, header

    if help or (not quiet) or (params and params.get_cc_mtz_pdb.verbose):
      print >>out,"Values of all params:"
      master_params.format(python_object=params).show(out=out)

    if help or params is None: return


    if (params.get_cc_mtz_pdb.pdb_in is None):
      print >>out,"\nSorry, please define an input model with pdb_in=mypdb.pdb"
      return
    self.pdb_path=params.get_cc_mtz_pdb.pdb_in
    if not os.path.isfile(params.get_cc_mtz_pdb.pdb_in):
       print >>out,"\nSorry, the file "+str(params.get_cc_mtz_pdb.pdb_in)+" is missing?"
       return
    if (params.get_cc_mtz_pdb.mtz_in is None):
       print >>out,"\nSorry you need an mtz file for get_cc_mtz_pdb"
       return
    if not os.path.isfile(params.get_cc_mtz_pdb.mtz_in):
       print >>out,"\nSorry the file "+str(params.get_cc_mtz_pdb.mtz_in)+" is missing?"
       return

    # do the work
    labin_1=self.set_up_labin_generic(
     mtz_file=params.get_cc_mtz_pdb.mtz_in,
     labin_input=params.get_cc_mtz_pdb.labin,
     require_phib=True,name='labin',target_map_coeffs=True)

    self.use_only_refl_present_in_mtz= \
       params.get_cc_mtz_pdb.use_only_refl_present_in_mtz

    self.verbose=params.get_cc_mtz_pdb.verbose


    self.scale=params.get_cc_mtz_pdb.scale
    if self.scale:
      print >>out,"Applying overall B factor and increment in B for each"
      print >>out,"atom beyond CB"

    self.Facts={'OutputDir':''}
    if params.get_cc_mtz_pdb.output_dir:
      self.Facts['OutputDir']=params.get_cc_mtz_pdb.output_dir

    if params.get_cc_mtz_pdb.temp_dir:
      if not os.path.exists(params.get_cc_mtz_pdb.temp_dir):
        os.mkdir(params.get_cc_mtz_pdb.temp_dir)
      self.Facts['temp_dir']=params.get_cc_mtz_pdb.temp_dir
    else:
      self.Facts['temp_dir']=self.create_temp_dir()

    self.remove_stopwizard()

    temp_pdb="TEMP_"+trim_file_name(params.get_cc_mtz_pdb.pdb_in).trimmed_file
    temp_mtz="TEMP_"+trim_file_name(params.get_cc_mtz_pdb.mtz_in).trimmed_file

    if params.get_cc_mtz_pdb.dry_run:
      print >>out, "ARGS: ",args
      return

    print >>out, "Get_cc_mtz_pdb: correlation of map and model allowing origin offsets"
    self.split_conformers=params.get_cc_mtz_pdb.split_conformers
    if self.split_conformers:
      print >>out,"Splitting conformers (separate analysis of each)"
      print >>out,"Note: requires all atoms in a conformer must be grouped together\n"
    self.fix_xyz=params.get_cc_mtz_pdb.fix_xyz
    self.rad_max=params.get_cc_mtz_pdb.rad_max
    if self.rad_max:
      print >>out,"CC radius (rad_max) will be fixed at %7.2f " %(self.rad_max)
      params.get_cc_mtz_pdb.fix_rad_max=True
    else:
      if params.get_cc_mtz_pdb.fix_rad_max:
        print >>out,"CC radius (rad_max) will be constant and set automatically"

    self.fix_rad_max=params.get_cc_mtz_pdb.fix_rad_max
    if self.fix_xyz:
       self.any_offset=False
       print >>out,"No offsets...fixing PDB coordinates"
    else:
      self.any_offset=params.get_cc_mtz_pdb.any_offset
      if self.any_offset:
         print >>out,"Allowing any offset -- even those not allowed by SG symmetry"

    copy_to_temp_dir(params.get_cc_mtz_pdb.temp_dir,
       params.get_cc_mtz_pdb.mtz_in,temp_mtz)

    if params.get_cc_mtz_pdb.atom_selection is None:
      copy_to_temp_dir(params.get_cc_mtz_pdb.temp_dir,
         params.get_cc_mtz_pdb.pdb_in,temp_pdb)
    else:  # apply atom selection and copy
      from phenix.autosol.get_pdb_inp import get_pdb_hierarchy
      hierarchy=get_pdb_hierarchy(file_name=params.get_cc_mtz_pdb.pdb_in)
      text=self.apply_atom_selection(
          params.get_cc_mtz_pdb.atom_selection,
          hierarchy=hierarchy,out=out,return_as_string=True)
      full_file=os.path.join(params.get_cc_mtz_pdb.temp_dir,temp_pdb)
      f=open(full_file,'w')
      print >>f, text
      f.close()
      print >>out,"Copied selection to %s " %(full_file)
    if params.get_cc_mtz_pdb.select_by_b:
      from phenix.autosol.get_pdb_inp import get_pdb_hierarchy
      hierarchy=get_pdb_hierarchy(file_name=params.get_cc_mtz_pdb.pdb_in)
      hierarchy = hierarchy.select(high_u_a)
      text=hierarchy.as_pdb_string()
      full_file=os.path.join(params.get_cc_mtz_pdb.temp_dir,temp_pdb)
      f=open(full_file,'w')
      print >>f, text
      f.close()
      print >>out,"Copied selection to %s " %(full_file)

    value=self.run_cc(
          mtz=temp_mtz,
          pdb=temp_pdb,
          labin=labin_1,
          resolution=params.get_cc_mtz_pdb.resolution,
          chain_type=params.get_cc_mtz_pdb.chain_type,
          quick=params.get_cc_mtz_pdb.quick,
          offset_pdb=params.get_cc_mtz_pdb.offset_pdb,
          out=out
           )
    if value: return

    print >>out, "overall CC: ",self.found_overall
    print >>out,"local CC: ",self.found_region
    from phenix.utilities import citations
    citations.show(source='get_cc_mtz_pdb',out=out)

  def special_cases(self,args):
    # special cases for input files so user doesn't need to specify:
    new_args=[]
    arg_use=None
    for arg in args:
      if (os.path.isfile(arg)):
        if arg[-3:]=='mtz':
          arg_use='get_cc_mtz_pdb.mtz_in='+arg
        elif arg[-3:]=='pdb':
          arg_use='get_cc_mtz_pdb.pdb_in='+arg
        else:
          arg_use=arg
      else:
        arg_use=arg
      if arg_use is not None:new_args.append(arg_use)
    return new_args

  def run_cc(self,mtz=None,pdb=None,labin=None,resolution=None,chain_type=None,
    quick=False,offset_pdb=None,out=sys.stdout):
    if not labin:
      labin="FP=FP PHIB=PHIM FOM=FOMM"
    print >>out, "Map from: ",mtz," using labin ",labin
    print >>out, "Model from: ",pdb
    self.mtz=mtz
    self.labin_line=labin
    self.quick=quick

    self.Facts['resolve_size']=""
    self.Facts['labin_FP_PHIB_FOM']=self.labin_line
    print >>out, "LABIN LINE: ",self.labin_line
    self.Facts['resolution']=resolution

    self.Facts['chain_type']=chain_type

    if self.Facts['resolution'] and self.Facts['resolution']<1000.:
       self.Facts['resolution_line']="resolution 1000. "+\
            str(self.Facts['resolution'])
    else:
       self.Facts['resolution_line']=""


    self.pdb=trim_file_name(pdb).trimmed_file

    if not self.fix_xyz:
      value=self.get_pdb_offset(out=out,offset_pdb=offset_pdb)
      if value: return value

    value=self.get_cc(out=out)
    if value: return value

    #delete_dir(self.Facts['temp_dir'],clean_up=True)
    #self.Facts['temp_dir']=""



  def get_pdb_offset(self,offset_pdb=None,out=sys.stdout):
    # offset our pdb file to match mtz
    print >>out, "Offsetting ",self.pdb," to match ",self.mtz
    print >>out, "Getting FC from self.pdb..."

    self.Facts['logfile']='fcalc.log'
    try:
      my_run_resolve=run_resolve(Facts=self.Facts,
       temp_dir=self.Facts['temp_dir'],
       size=self.Facts['resolve_size'],
       hklin=trim_file_name(self.mtz).trimmed_file,
       labin=self.Facts['labin_FP_PHIB_FOM'],
       resolution_line=self.Facts['resolution_line'],
       command_3='fcalc '+str(self.pdb),
       logfile=self.Facts['logfile'])
    except KeyboardInterrupt: raise
    except:
      print >>out, "\nSorry, unable to get FC..."
      print >>out, "Please check that",self.Facts['labin_FP_PHIB_FOM'],\
        "specifies your column labels correctly"
      print >>out, "Please also check that",self.pdb_path," is ok"
      return 1

    self.fcalc=os.path.join(self.Facts['temp_dir'],'resolve.mtz')
    if not os.path.isfile(self.fcalc):
      raise Sorry( "Sorry, unable to get FC from %s " %(self.mtz))

    print >>out, "FC is in ",self.fcalc

    # now offset self.pdb to match self.mtz:

    self.pdb_out=os.path.join(self.Facts['temp_dir'],'resolve.pdb')
    if os.path.isfile(self.pdb_out): os.remove(self.pdb_out)

    self.Facts['logfile']='offset.log'
    if self.any_offset:
      command_2="any_offset"
    else:
      command_2=""
    try:
      my_run_resolve=run_resolve(Facts=self.Facts,
       temp_dir=self.Facts['temp_dir'],
       size=self.Facts['resolve_size'],
       hklin=trim_file_name(self.mtz).trimmed_file,
       labin=self.Facts['labin_FP_PHIB_FOM'],
       hklperfect=trim_file_name(self.fcalc).trimmed_file,
       labperfect="FP=FP PHIB=PHIM FOM=FOMM",
       resolution_line=self.Facts['resolution_line'],
       command_1='hkl_offset_file offset.mtz',
       command_2=command_2,
       command_3='pdb_in '+str(self.pdb),
       logfile=self.Facts['logfile'])

    except KeyboardInterrupt: raise
    except Exception,e:
      raise Sorry("\nSorry, failed... \n %s" %(str(e)))

    if not os.path.isfile(self.pdb_out):
      raise Sorry( "Sorry, unable to get offset pdb from %s" %(self.pdb,self.mtz))

    self.pdb_out_new=offset_pdb

    from phenix.autosol.copy_extra_to_pdb import copy_extra_to_pdb
    from libtbx.utils import null_out
    copy=copy_extra_to_pdb(
        pdb1=self.pdb_out,
        pdb2=os.path.join(self.Facts['temp_dir'],self.pdb),
        pdb3=os.path.join(self.Facts['temp_dir'],self.pdb_out_new),
        copy_coords=True,out=null_out())
    copy.run()

    print >>out, "Offset pdb file is in ",self.pdb_out_new
    copy_from_temp_dir(self.Facts['temp_dir'],self.pdb_out_new,self.pdb_out_new,
         OutputDir=self.Facts['OutputDir'])
    copy_from_temp_dir(self.Facts['temp_dir'],offset_pdb,
              offset_pdb,OutputDir=self.Facts['OutputDir'])



  def get_cc(self,out=sys.stdout):  # offset our pdb file to match mtz
    if self.fix_xyz:
      print >>out, "Getting CC of ",self.pdb," relative to ",self.mtz
      self.pdb_out_new=self.pdb
    else:
      print >>out, "Getting CC of ",self.pdb," (offset to ",self.pdb_out_new,\
         ") to match ",self.mtz

    self.Facts['logfile']='cc.log'
    resolution=self.Facts['resolution']
    delete_file(self.Facts['logfile'])
    delete_file(os.path.join(self.Facts['temp_dir'],self.Facts['logfile']))
    if self.scale:
      command_1='evaluate_model'
    else:
      command_1='evaluate_model \nb_overall 0'
    if not self.use_only_refl_present_in_mtz:
      command_1+="\nfill_ratio 1.0 \nfill \nres_fill "+str(resolution)

    if self.split_conformers:
      command_1+="\nconformers"
    if self.rad_max:
      command_1+="\nrad_max %7.2f" %(self.rad_max)
    if self.fix_rad_max:
      command_1+="\nfix_rad_max"
    if self.verbose:
     command_1+="\nverbose"
    if self.quick: command_1+="\nno_cc"
    my_run_resolve=run_resolve(Facts=self.Facts,
       temp_dir=self.Facts['temp_dir'],
       size=self.Facts['resolve_size'],
       hklin=trim_file_name(self.mtz).trimmed_file,
       labin=self.Facts['labin_FP_PHIB_FOM'],
       resolution_line=self.Facts['resolution_line'],
       command_1=command_1,
       command_3='model '+str(trim_file_name(self.pdb_out_new).trimmed_file),
       logfile=self.Facts['logfile'])
    copy_from_temp_dir(self.Facts['temp_dir'],self.Facts['logfile'],
                  self.Facts['logfile'],OutputDir=self.Facts['OutputDir'])

    print >>out, "Detailed analysis of correlation of map and model are in: ",\
        self.Facts['logfile']
    if not self.fix_xyz:
      print >>out, "Offset PDB file is in %s\n" %(self.pdb_out_new)
    self.found_overall=0
    self.found_region=0
    for line in open(os.path.join(self.Facts['temp_dir'],self.Facts['logfile'])).readlines():
       if not self.found_overall and re.search('Overall map correlation:',line):
         try:
           self.found_overall=string.atof(string.split(line)[3])
         except KeyboardInterrupt: raise
         except: pass
       if not self.found_region and \
         re.search('Map correlation in region of model:',line):
         try:
           self.found_region=string.atof(string.split(line)[6])
         except KeyboardInterrupt: raise
         except: pass
    print >>out, "Correlation in region of model: ",self.found_region,\
         "...overall: ",self.found_overall


  def create_temp_dir(self):
      import os
      # make a directory that does not exist TEMPxx
      for i in xrange(1000):
        temp_dir='TEMP'+str(i)
        if not os.path.exists(temp_dir):
          os.mkdir(temp_dir)
          return temp_dir

      line="Failed to make TEMP directory..."
      raise AssertionError,line


  def get_summary_and_header(self,command_name):
    header="\n"
    header+="\n#                       "+str(command_name)
    header+="\n#"
    header+="\n# Get correlation between atoms in a PDB file and map "
    header+="\n# offsetting the PDB file by allowed origin shifts"
    header+="\n\n# Type phenix.doc for help"

    summary= "usage: phenix.%s protein.pdb mtzfile.mtz [labin='FP=FP PHIB=PHIM FOM=FOMM'] " % command_name
    return summary,header

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    os.makedirs(self.output_dir)
    os.chdir(self.output_dir)
    return get_cc_mtz_pdb(args=list(self.args), out=sys.stdout)

def validate_params (params, callback=None) :
  if None in [params.get_cc_mtz_pdb.mtz_in,
              params.get_cc_mtz_pdb.labin] :
    raise Sorry("Missing the map file and/or data labels.")
  elif params.get_cc_mtz_pdb.pdb_in is None :
    raise Sorry("Please supply a PDB file.")

def get_cc_values (results) :
  (overall_cc, local_cc) = (None, None)
  if hasattr(results, "found_overall") :
    overall_cc = results.found_overall
  if hasattr(results, "found_region") :
    local_cc = results.found_region
  return (overall_cc, local_cc)

def finish_job (results) :
  stats = []
  (overall_cc, local_cc) = get_cc_values(results)
  stats.append(("Overall CC", format_value("%.4f", overall_cc)))
  stats.append(("Local CC", format_value("%.4f", local_cc)))
  return ([], stats)

if (__name__ == "__main__"):
  argument_list=sys.argv[1:]
  if True or is_debug(argument_list).value:

    sys_stdout_sav=sys.stdout
    get_cc_mtz_pdb=get_cc_mtz_pdb(args=sys.argv[1:])
  else:
   try:
    sys_stdout_sav=sys.stdout
    get_cc_mtz_pdb=get_cc_mtz_pdb(args=sys.argv[1:])
   except KeyboardInterrupt:
    pass
   except Exception, e:
    print "\n************************************************"
    print e
    print "\n************************************************"
    if sys.stdout != sys_stdout_sav:
      sys.stdout=sys_stdout_sav # restore output stream so we can print error
      print  "\n************************************************"
      print e
      print  "\n************************************************"
    from phenix.utilities.is_raise_sorry import is_raise_sorry
    if is_raise_sorry(argument_list).value:
      from libtbx.utils import Sorry
      raise Sorry(e)

