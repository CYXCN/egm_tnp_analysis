import copy

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

def createBins(bining, cut, to_correct_cut=None):
    """
    Create bins from bining definition.
    Supports distinguishing between standard variables and corrected variables.
    
    bining: list of dicts or dict of lists (regions)
    cut: Base cut string for standard variables (will be extracted as baseSelection)
    to_correct_cut: Base cut string for variables that need correction
    """

    # --- 1. Handle multiple regions (Dictionary Mode) ---
    if isinstance(bining, dict):
        all_bins = []
        for region_name, region_bining in bining.items():
            # Recursive call, passing to_correct_cut
            region_bindef = createBins(region_bining, cut, to_correct_cut)
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

    # --- 2. Handle single list (List Mode) ---
    
    # Initialize N-dimensional binning indices
    nbin = 1
    index = list(range(len(bining)))
    for ix in range(len(index)):
        index[ix] = -1
    listOfIndex = []    
    listOfIndex.append(index)

    ### Preprocessing: calculate total number of bins and generate index combinations
    for iv in range(len(bining)):
        # Identify variable name (var or cor_var)
        if 'cor_var' in bining[iv]:
            var_raw = bining[iv]['cor_var']
        else:
            var_raw = bining[iv]['var']

        # Check configuration completeness
        if 'type' not in bining[iv] or 'bins' not in bining[iv]:
            print(('bining is not complete for var %s' % var_raw))
            # Fallback return with base cuts
            return {
                'vars': [],
                'bins': [{'cut': None, 'cor_cut': parse_and_to_list(to_correct_cut) if to_correct_cut else []}],
                'baseSelection': cut
            }

        bins_spec = bining[iv]['bins']
        nb1D = 1
        
        # Detect pair-mode
        if len(bins_spec) > 0 and isinstance(bins_spec[0], (list, tuple)):
            nb1D = len(bins_spec)
        else:
            if bining[iv]['type'] == 'float':
                nb1D = len(bins_spec) - 1
            elif bining[iv]['type'] == 'int':
                nb1D = len(bins_spec)

        nbin = nbin * nb1D

        # Expand index list (Cartesian product logic)
        listOfIndexInit = copy.deepcopy(listOfIndex)
        for ib_v in range(nb1D):
            if ib_v == 0:
                for ib in range(len(listOfIndex)):
                    listOfIndex[ib][iv] = ib_v            
            else: 
                for ib in range(len(listOfIndexInit)):
                    listOfIndexInit[ib][iv] = ib_v
                listOfIndex.extend(copy.deepcopy(listOfIndexInit))

    listOfBins = []
    ibin = 0
    nbins = len(listOfIndex)

    ### Generate specific bin definitions
    for ix in listOfIndex:
        # Initialize current bin variables - no longer include base cut here
        binCut = None  # Only store bin-specific cuts
        binCorCut = None  # New: for corrected variable cuts
        
        # Do NOT add base cut here - it will be in baseSelection
        if to_correct_cut is not None:
            binCorCut = to_correct_cut

        # Generate bin name prefix
        if nbins > 10000: binName = 'bin%d' % ibin
        elif nbins > 1000: binName = 'bin%04d' % ibin
        elif nbins > 100: binName = 'bin%03d' % ibin
        else: binName = 'bin%02d' % ibin

        binTitle = ''
        binVars = {}

        # Iterate over each dimension (variable) of this bin
        for iv in range(len(ix)):
            # 1. Identify variable type and name
            is_correction_var = False
            if 'cor_var' in bining[iv]:
                var_raw = bining[iv]['cor_var']
                is_correction_var = True
            else:
                var_raw = bining[iv]['var']
                is_correction_var = False

            # Handle abs()
            if isinstance(var_raw, str) and var_raw.startswith('abs(') and var_raw.endswith(')'):
                var_name = var_raw[4:-1]
                var_expr = f'abs({var_name})'
                var_key = f'abs_{var_name}'
            else:
                var_name = var_raw
                var_expr = var_raw
                var_key = var_raw

            bins1D = bining[iv]['bins']
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
            
            # --- Int type handling ---
            elif varType == 'int':
                val = bins1D[ix[iv]]
                minv, maxv = val, val
                cond = f'{var_expr} == {val}'
                title_part = f'{var_expr} = {val}'
                name_part = f'{var_key}Eq{val}'

            # --- Assemble cut string ---
            if cond:
                # Key modification: route to different cut strings based on whether it's a cor_var
                if is_correction_var:
                    if binCorCut is None:
                        binCorCut = cond
                    else:
                        binCorCut = f'{binCorCut} && {cond}'
                else:
                    if binCut is None:
                        binCut = cond
                    else:
                        binCut = f'{binCut} && {cond}'
                
                # Title and Name should always be updated
                if binTitle == '':
                    binTitle = title_part
                else:
                    binTitle = f'{binTitle}; {title_part}'
                
                binName = f'{binName}_{name_part}'

            binVars[var_key] = {'min': minv, 'max': maxv}

        # Clean up name string
        binName = binName.replace('-', 'm')
        binName = binName.replace('.', 'p')

        # Append to result list, adding cor_cut field
        listOfBins.append({
            'cut': binCut,  # Now only contains bin-specific cuts
            'cor_cut': parse_and_to_list(binCorCut), 
            'title': binTitle, 
            'name': binName, 
            'vars': binVars
        })
        ibin += 1

    # --- Collect variable list (for 'vars' in return value) ---
    listOfVars = []
    for iv in range(len(bining)):
        if 'cor_var' in bining[iv]:
            var_raw = bining[iv]['cor_var']
        else:
            var_raw = bining[iv]['var']
            
        if isinstance(var_raw, str) and var_raw.startswith('abs(') and var_raw.endswith(')'):
            var_name = var_raw[4:-1]
            listOfVars.append(f'abs_{var_name}')
        else:
            listOfVars.append(var_raw)
    
    binDefinition = {
        'vars': listOfVars,
        'bins': listOfBins,
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
    

