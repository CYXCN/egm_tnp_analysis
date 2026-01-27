import copy

# help function to parse cut strings into list of variable conditions
def parse_and_to_list(expr_str):
    result = []
    conditions = expr_str.split(' && ')
    
    for cond in conditions:
        var_name = ""
        min_val = -9999
        max_val = 9999
        
        if '>=' in cond:
            var_name, val = cond.split('>=')
            min_val = float(val)
        elif '>' in cond:
            var_name, val = cond.split('>')
            min_val = float(val) + 1e-9
        elif '<=' in cond:
            var_name, val = cond.split('<=')
            max_val = float(val)
        elif '<' in cond:
            var_name, val = cond.split('<')
            max_val = float(val) - 1e-9
        else:
            var_name = cond.strip()
        
        var_name = var_name.strip()
        existing_item = next((item for item in result if item['var_name'] == var_name), None)
        
        if existing_item:
            existing_item['min'] = max(existing_item['min'], min_val)
            existing_item['max'] = min(existing_item['max'], max_val)
        else:
            result.append({
                'var_name': var_name,
                'min': min_val,
                'max': max_val
            })
    
    return result

def createBins(bining, cut, to_correct_cut=None, legacy_mode=False):
    """
    Create bins from bining definition.
    Supports distinguishing between standard variables and corrected variables.
    
    bining: list of dicts or dict of lists (regions)
    cut: Base cut string for standard variables (will be extracted as baseSelection)
    where the cor_var means the variable that needs correction
    to_correct_cut: Base cut string for variables that need correction
    legacy_mode: If True, includes baseSelection in each bin's cut (for backward compatibility)
    """

    # --- 1. Handle multiple regions (Dictionary Mode) ---
    # iterate over each region and combine results
    if isinstance(bining, dict):
        all_bins = []
        for region_name, region_bining in bining.items():
            # Recursive call, passing to_correct_cut and legacy_mode
            region_bindef = createBins(region_bining, cut, to_correct_cut, legacy_mode=legacy_mode)
            # Add region name to bin names
            for bin_info in region_bindef['bins']:
                bin_info['region'] = region_name
            all_bins.extend(region_bindef['bins'])
        
        # Renumber all bins with unified naming
        total_bins = len(all_bins)
        for idx, bin_info in enumerate(all_bins):
            # Keep original naming logic, preserve original suffix
            info_part = bin_info['name'].split('_', 1)[-1]  
            if total_bins > 100:
                prefix = f'bin{idx:03d}'
            elif total_bins > 1000:
                prefix = f'bin{idx:04d}'
            elif total_bins > 10000:
                prefix = f'bin{idx}'
            else:
                prefix = f'bin{idx:02d}'
            bin_info['name'] = f'{prefix}_{info_part}'
        
        # Collect all variable names (including var and cor_var)
        all_vars = []
        for region_bining in bining.values():
            for var_def in region_bining:
                # Prefer checking cor_var first
                if 'cor_var' in var_def:
                    var_raw = var_def['cor_var']
                else:
                    var_raw = var_def['var']
                
                if isinstance(var_raw, str) and var_raw.startswith('abs(') and var_raw.endswith(')'):
                    var_name = var_raw[4:-1]
                    var_key = f'abs_{var_name}'
                else:
                    var_key = var_raw
                if var_key not in all_vars:
                    all_vars.append(var_key)
        
        return {
            'vars': all_vars,
            'bins': all_bins,
            'baseSelection': cut  # Add baseSelection at top level
        }

    # --- 2. Handle single list, original mode ---
    nbin = 1
    index = list(range(len(bining)))
    for ix in range(len(index)):
        index[ix] = -1
    listOfIndex = []    
    listOfIndex.append( index )

    ### first map nD bins in a single list
    for iv in range(len(bining)):
        is_correction_var = True if 'cor_var' in bining[iv] else False
        var = bining[iv]['var'] if not is_correction_var else bining[iv]['cor_var']
        if 'type' not in bining[iv] or 'bins' not in bining[iv]:
            print(('bining is not complete for var %s' % var))
            return listOfIndex
        nb1D = 1
        if   bining[iv]['type'] == 'float' :
            nb1D = len(bining[iv]['bins'])-1
        elif bining[iv]['type'] == 'int' :
            nb1D = len(bining[iv]['bins'])
        # Detect pair-mode, where bins are defined as pairs of edges
        if len(bining[iv]['bins']) > 0 and isinstance(bining[iv]['bins'][0], (list, tuple)):
            nb1D = len(bining[iv]['bins'])
        nbin = nbin * nb1D

        listOfIndexInit = copy.deepcopy(listOfIndex)
        for ib_v in range(nb1D):
            if ib_v == 0 :
                for ib in range(len(listOfIndex)):
                    listOfIndex[ib][iv] = ib_v            
            else: 
                for ib in range(len(listOfIndexInit)):
                    listOfIndexInit[ib][iv] = ib_v
                           
                listOfIndex.extend(copy.deepcopy(listOfIndexInit))

    listOfBins = []
    ibin = 0
    nbins = len(listOfIndex)
    for ix in listOfIndex:
        ### make bin definition
        binCut   = None
        binCorCut = to_correct_cut  if to_correct_cut is not None else None
        binName  = 'bin%02d'%ibin
        if nbins > 100   :  binName  = 'bin%03d'%ibin
        if nbins > 1000  :  binName  = 'bin%04d'%ibin
        if nbins > 10000 :  binName  = 'bin%d'%ibin

        binTitle = ''
        binVars = {}

        for iv in range(len(ix)):
            is_correction_var = True if 'cor_var' in bining[iv] else False
            var     = bining[iv]['var'] if not is_correction_var else bining[iv]['cor_var']

            # Handle abs()
            if isinstance(var, str) and var.startswith('abs(') and var.endswith(')'):
                var_name = var[4:-1]
                var_expr = f'abs({var_name})'
                var_key = f'abs_{var_name}'
            else:
                var_name = var
                var_expr = var
                var_key = var

            bins1D  = bining[iv]['bins']
            varType = bining[iv]['type']
            
            # pair-mode detection
            pairMode = (len(bins1D) > 0 and isinstance(bins1D[0], (list, tuple)))
            
            cond = None      # Condition string for current dimension
            title_part = ''  # Title part for current dimension
            name_part = ''   # Name part for current dimension
            minv, maxv = None, None

            # --- Float type handling ---
            if varType == 'float':
                if pairMode:
                    minv = bins1D[ix[iv]][0]
                    maxv = bins1D[ix[iv]][1]
                    
                    if minv is None and maxv is None:
                        cond = None
                    elif minv is None:
                        cond = f'{var_expr} < {maxv}'
                        title_part = f'{var_expr} < {maxv:.3f}'
                        name_part = f'{var_key}_To{maxv:.3f}'
                    elif maxv is None:
                        cond = f'{var_expr} >= {minv}'
                        title_part = f'{var_expr} >= {minv:.3f}'
                        name_part = f'{var_key}_{minv:.3f}ToInf'
                    else:
                        cond = f'{var_expr} >= {minv} && {var_expr} < {maxv}'
                        title_part = f'{minv:.3f} < {var_expr} < {maxv:.3f}'
                        name_part = f'{var_key}_{minv:.3f}To{maxv:.3f}'
                else:
                    # Original edge-based behavior
                    minv = bins1D[ix[iv]]
                    maxv = bins1D[ix[iv]+1]
                    cond = f'{var_expr} >= {minv} && {var_expr} < {maxv}'
                    title_part = f'{minv:.3f} < {var_expr} < {maxv:.3f}'
                    name_part = f'{var_key}_{minv:.3f}To{maxv:.3f}'

                
            if varType == 'int' :
                val = bins1D[ix[iv]]
                minv, maxv = val, val
                cond = f'{var_expr} == {val}'
                title_part = f'{var_expr} = {val}'
                name_part = f'{var_key}Eq{val}'

            # --- Assemble cut string ---
            if cond:
                # Key modification: route to different cut strings based on whether it's a cor_var
                if is_correction_var:
                    binCorCut = cond if binCorCut is None else f'{binCorCut} && {cond}'
                else:
                    binCut = cond if binCut is None else f'{binCut} && {cond}'
                binTitle = title_part if binTitle == '' else f'{binTitle}; {title_part}'
                binName = f'{binName}_{name_part}'

            binVars[var_key] = {'min': minv, 'max': maxv}

            binName = binName.replace('-','m')
            binName = binName.replace('.','p')

        listOfBins.append({
            'cut': f'{cut} && {binCut}' if legacy_mode and cut and binCut else binCut,  # Legacy mode: include baseSelection
            'cor_cut': parse_and_to_list(binCorCut), 
            'title': binTitle, 
            'name': binName, 
            'vars': binVars
        })
        ibin = ibin + 1
 
    listOfVars = []
    for iv in  range(len(bining)):
        if 'cor_var' in bining[iv]:
            var = bining[iv]['cor_var']
        else:
            var = bining[iv]['var']            
        if isinstance(var, str) and var.startswith('abs(') and var.endswith(')'):
            var = 'abs_' + var[4:-1]
        listOfVars.append(var)
        
    binDefinition = {
        'vars' : listOfVars,
        'bins' : listOfBins,
        'baseSelection': cut  # Add baseSelection as independent field
        } 
    return binDefinition

def tuneCuts( bindef, cuts ) :
    if cuts is None:
        return
    
    for ibin in list(cuts.keys()):
        cut0 = bindef['bins'][ibin]['cut']
        cut1 = cuts[ibin]
        bindef['bins'][ibin]['cut'] = '%s && %s ' % (cut0,cut1)
    

