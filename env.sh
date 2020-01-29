module unload specsim
module unload desisim


export  SPECSIM_ROOT=/global/homes/m/mjwilson/desi/bluebump/

export  DESIMODEL=/global/homes/m/mjwilson/desi/bluebump/

export  PYTHONPATH=$SPECSIM_ROOT/specsim:$SPECSIM_ROOT/specsim/py:$PYTHONPATH                                                                                                          
export  PATH=$SPECSIM_ROOT/specsim/bin:$PATH  

export  PYTHONPATH=$SPECSIM_ROOT/desisim:$SPECSIM_ROOT/desisim/py:$SPECSIM_ROOT/desisim/scripts/:$PYTHONPATH
export  PATH=$SPECSIM_ROOT/desisim/bin:$PATH
