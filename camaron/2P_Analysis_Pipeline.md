# 2P Analysis Pipeline

*Github Repository → ImagingCode-Glickfeld-Hull\camaron*

1. **oriAdapt_V1_cam.m** # structure containing info for each experiment
2. **SingleChannelTC_2AFC_OriAdapt_udpated.m** AND **SingleChannelTC_2AFC_OriAdapt_passive_only.m** (for naive mice) 
    1. Doing many things...see below:
    Input -> raw imaging and behavior data for a single session
    Output -> Behaving, Direction Tuning, Passive info
    2. Used in conjunction with **RegistrationPathways.m (PUSH TO GIT)** for days with signal interference

```markdown
 # Behaving Image
 _reg_shifts.mat ('out', 'data_avg', 'diff_out', 'move_ind')
 _input.mat ('adapt_input')
 _redData.mat ('green_data_avg', 'red_data_avg')
 _photoData.mat ('photoLoc', 'stimOnFrames', 'cStimOn','cAdaptOn') OR ('stimOnFrames', 'cStimOn','cAdaptOn') 
**Depends on photoFrameFinder (ephys OR Sandworks)**

 # Dir Tuning
 _reg_shifts.mat ('stimOnFrames', 'cStimOn','cAdaptOn')
 _input.mat ('dir_input')
 _stimData.mat ( 'dir_mat', 'Dirs', 'nDir', 'dir_dfof_all')

 # Cell Segmentation & Neuropiul Subtraction
 _mask_cell.mat ('data_dfof_max', 'mask_cell')
 _TCs.mat ('data_tc','np_tc','npSub_tc')

 # Adapt Analysis (Behavior)
 _adaptResp.mat ('mask_label', 'data_adapt_dfof', 'adapt_cyc_resp', 'base_win_all', 'resp_win_all','adapt_resp_ind')
 _stimResp.mat ('stim_resp', 'stim_resp_ind', 'base_win', 'resp_win) 
			# data_stim_dfof not saved!!!!
     > data_stim_dfof = bsxfun(@rdivide, bsxfun(@minus, data_stim, dataf), dataf);
     > Generated at line 658 (Create loop to append)
		 > Append with **append_data_stim_dfof.m** (intial TC script corrected!!!)

 _stim_Data.mat ('tGratingOri', 'tOris', 'aGratingOri', 'aOris', 'aGratingContrast', 'ind_cond', 'SIx', 'MIx')

 # Dir Tuning NP Subtraction and Analysis
 _TCs.mat ('data_tc', 'np_tc', 'npSub_tc_dir')
 _oriTuningandFits.mat ('avgResponseEaOri','semResponseEaOri','vonMisesFitAllCellsAllBoots','fitReliability','R_square', 'tuningTC')
 -oriTuningInfo.mat ('prefOri', 'prefOri_bootdiff', 'ind_theta90', 'tunedCells', 'edges')

 # Passive Image, NP Subtraction
 _reg_shifts.mat ('out', 'data_avg')
 _input.mat ('pass_input')
 _TCs.mat ('data_tc','np_tc','npSub_tc)
 _photoData.mat ('photoLoc', 'stimOnFrames', 'cStimOn','cAdaptOn') OR ('stimOnFrames', 'cStimOn','cAdaptOn')

 # Adapt Analysis (Passive)
 _adaptResp.mat ('mask_label', 'data_adapt_dfof', 'adapt_cyc_resp', 'base_win_all', 'resp_win_all','adapt_resp_ind')
 _stimResp.mat ('stim_resp',  'data_stim_dfof', 'stim_resp_ind', 'base_win', 'resp_win)
 _stim_Data.mat ('tGratingOri', 'tOris', 'aGratingOri', 'aOris', 'aGratingContrast', 'ind_cond', 'SIx')

# Many unlisted, intermediate figures (.pdf)

```

1. **genExptListFromStruct.m** # Get good TC extractions with including behaving, passive, direction tuning, FB OFF, and s > 0.8
    1. **PUSH TO GIT**
    2. Input - > Data struct, behavior input.mat
    3. Output → good_expt_list.mat (mouse_expts, mouse_good_list, …, …)
2. **FOV_analysis.m**  # Plot image FOVs to view imaging sessions at different depths with non-overlapping and responsive cells for each mouse
    1. **PUSH TO GIT**
    2. Input → good_expt_list.mat, raw image files, _redData.mat, _mask_cell.mat, adaptResp.mat, stimResp.mat
        1. intermediate → adapt_and_stim_resp_ind = union(adapt_resp_ind{:}, stim_resp_ind);
    3. output → 
        1. mouse FOV_analysis.pdf
        2. mouse _image_matching.mat (expt_list_final)
3. **pool_2AFC_2P_data.m** # pool 2AFC responses by mouse and create structure containing everything needed for final analysis 
    1. Inputs → TC_extraction outputs + expt_list_final
    2. Outpus → Single data structure:
        1. mouse_pooled_2AFC_2P_data.mat (expt_list_final, data_adapt_dfof_all, data_stim_df_all, adapt_and_stim_resp_ind_all, ori_bins …. …. …. … ) 
4. **figs_2AFC_2P.m** # Post processing data and plotting figures for a single mouse…
    1. Inputs → mouse_pooled_2AFC_2P_data.mat
    2. Figures….