#!/usr/bin/env python3

# %% ========================================================================
# % 
# %                          peakFitting function
# %  Yoon Research Group, Deparment of Geological Science, University of Florida
# % -------------------------------------------------------------------------
# % This code is to automate the multi-peak fitting process needed for 
# % anlayzing the curves produced by a spectrofluorophotometer. The code 
# % estimates the peak area correponding to a dye tracer. 
# % 
# % Currently, this code can handle uranine (fluorescein), sulforhodamine B 
# % (SrB), and rhodamine WT (RWT) which typically have the ranges of peak 
# % emission wavelength range as follows:
# %   - when dissolved in watter
# %       * uranine: 506 - 510 nm
# %       * eosin: 532-538 nm
# %       * RWT: 572 - 577 nm 
# %       * SrB: 580 - 583 nm 
# % 
# % This peak fitting code is based the PYTHON using non-linear optimization 
# % function 'optimize', method "SLSQP"
# % with the contraint variabels A, B, Aeq, Beq, lb,& ub. 
# % dye_type = currently only 'fluorescein', "SrB" 
# % sample_type =  either 'water' or 'eluant'
# % =========================================================================


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import os

# =========================================================================
# Sub-functions
# =========================================================================

def pearson7(x, params):
    """
    Pearson VII function.
    Parameters:x (array): Input data, (wavelength)
        params (array): Parameters of the Pearson VII function.
    Returns: (array) Output of the Pearson VII function, (intensity simulated)
    """
    return params[1] / (1 + ((x - params[0]) / params[2]) ** 2 * (2 ** (1 / params[3]) - 1)) ** params[3]

def peak_sum(x, params):
    """
    Sum of multiple Pearson VII functions.
    Parameters:x (array): Input data (wavelength)
        params (array): Parameters of the Pearson VII functions (center, high, hwhm,shape)
    Returns: (array): Output of the sum of Pearson VII functions.
    """
    n_peak = len(params) // 4
    curve = pearson7(x, params[:4])
    for i in range(1, n_peak):
        curve= curve + pearson7(x, params[i * 4:(i + 1) * 4])
    return curve

# =========================================================================
# Main-function
# =========================================================================

def peak_fitter_so(data_dir, fig_dir, fl_name, n_iter):
    """
    Peak fitter function.
    Parameters:
        data_dir (str): Directory of the data file.
        fig_dir (str): Directory of the figure file.
        file_name (str): Name of the data file.
        sample_type (str): Type of the sample (water,eluent)
        n_iter (int): Number of iterations.
    Returns:
        area (list): The list of calculated peak areas.
    """
    file_path = os.path.join(data_dir, fl_name)  ### Construct the full path to the file
    wavelength = []      # Initialize the variables
    intensity = []

    try:
        # Open the file and read its contents
        with open(file_path, 'r') as fid:
            lines = fid.readlines()

        # Find the data start line (search for 'Wavelength nm')
        data_start = 0
        for iL, line in enumerate(lines):
            if 'Wavelength nm' in line:
                data_start = iL + 1  # Start reading data after this line
                break

        # Process the data lines after the header
        for line in lines[data_start:]:
            note = line.split(',')  ### Split the line by comma and extract wavelength and intensity values
            wavelength.append(float(note[0]))
            intensity.append(float(note[1]))

    except FileNotFoundError:
        print(f"Error: The file at {file_path} was not found.")
        return []
    except ValueError as e:
        print(f"Error while reading the file: {e}")
        return []

    # Convert lists to numpy arrays for easier mathematical processing
    wavelength = np.array(wavelength)
    intensity = np.array(intensity)
    
    def cost_fn(params):
        return np.sqrt(np.mean((intensity - peak_sum(wavelength, params)) ** 2))
    
    
    ##Matrix with NOMS and dyes          
    ##center  uranine (505-509), Eosin (532-537), RWT (571-577), SrB (580 - 583) (POSITION 8,9,10,11) from index-=0
    params = np.array([
    np.linspace(300, 500, 7).tolist() + [540]+ [508, 534, 572, 582] + [600, 640, 680, 730, 760],  # center
    [10] * 8 + [500, 500, 500, 500] + [10] * 5,  # height
    [100] + [30] * 7 + [10, 10,10, 10] + [30] * 5,  # hwhm
    [2] * 8 + [10, 10,10, 10] + [30] * 5  # shape
    ])   
    params=params.ravel(order='F') ## order [c1,h1,hw1,s1, c2,h2,hw2,s2,,,,c17,h17,hw17,s17]
    n_peak = len(params) // 4   ##18 #natural organic materials

    fixed = np.array([
        [0] * n_peak,  # center
        [0] * n_peak,  # height
        [0] * n_peak,  # hwhm
        [0] + [1] * (n_peak - 1)  # shape
    ])
    fixed=fixed.ravel(order='F')
    
    shape = np.array([
        [0] * n_peak,  # center
        [0] * n_peak,  # height
        [0] * n_peak,  # hwhm
        [0] + [2] * 7 + [10, 10,10, 10] + [5] * 5  # shape
    ])        
    shape=shape.ravel(order='F')
            
    Aeq = np.diag(fixed) ### it should diag the vector (fixed) from the estructure [c1,h1,hwhm,s1...c2,h2,hw2,s3.... s17] (68x68)
    beq = shape
    
    lb = np.array([
        [320, 340, 360, 370, 410, 430, 480, 525]+[506, 532, 555, 580] + [600] * 5,  # center
        [0] * n_peak,  # height
        [5] * n_peak,  # hwhm
        [0.4] + [0]*7+ [8.5,8.5,8.5,8.5]+ [0]*(n_peak - 12)  # shape
    ])
    lb=lb.ravel(order='F') 

    ub = np.array([
        [380, 400, 450, 430, 470, 490, 500, 545]+ [510, 538, 565, 584] + [850] * 5,  # center
        [600] * 7 + [350]+ [np.inf, np.inf, np.inf, np.inf] + [500] * 5,  # height
        [300] * 8 + [15,15,15,15] + [100] * 5,  # hwhm
        [0.9] + [5]*7+ [10.1,10.1,10.1,10.1] +[100] *5  # shape
    ])
    ub=ub.ravel(order='F') 

        
##### optimization processe
    op_p=[]  # To store optimized params for each series  
    for i_iter in range(n_iter):  #iter =0 
        print(f'{i_iter+1}/{n_iter} iteration')
        
        if i_iter >0:  # Use the previous optimized parameters as the initial guess for the next iteration
            params = op_p[-1]  # Use the last set of optimized parameters
        else:
            if len(op_p) > 0:   # For the first iteration, check if op_p has any previous values
                params = op_p[-1]  # Use the last optimized parameters from previous dataset
            else:
                params = params  # Use your initial guess for the first dataset iteration 
        
        def eq_constraint(params):
            return np.dot(Aeq, params) - beq
        
        def filter_infs(lb, ub):
            return [(l if np.isfinite(l) else -1e6, u if np.isfinite(u) else 1e6) for l, u in zip(lb, ub)]
        bounds = filter_infs(lb, ub)
        
        constraints = [{'type': 'ineq', 'fun': eq_constraint}]  ##ineq
        res=minimize(cost_fn, params, method='SLSQP', bounds=bounds, constraints=constraints,
                      options = {'disp': False})
        params =res.x ##update last parameters optimized
        op_p.append(params)  # Store the optimized parameters for later use
        curve = peak_sum(wavelength, params)
        
        # Visualization
        plt.figure(figsize=(8.5,8))

        # Plotting
        from matplotlib.gridspec import GridSpec
        gs = GridSpec(2, 1, height_ratios=[1.5,0.8])  # Give row 1 (second subplot) 3x the height of the others
        
        #goodness of fit metrics
        SSR = np.sum((intensity-curve) ** 2) # sum of squeared residual
        R2 = 1 - SSR / np.sum((intensity - np.mean(intensity)) ** 2) #Coef of determinantion

    
        # === Top Subplot ===
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        ax1 = plt.subplot(gs[0])
        ax1.plot(wavelength, intensity, label='Measured', color="black", linestyle=':', linewidth=3)
        ax1.plot(wavelength, curve, label='Fitted', color=colors[1] )
        ax1.legend(fontsize=16)
        ax1.set_ylabel('Intensity [-]', fontsize=16)
        ax1.set_xticklabels([])  #Remove x-axis numbers
        ax1.set_xlabel(None)  # Remove xlabel if present
        ax1.tick_params(labelsize=16, width=2)  # Increase font size of axis numbers

        ###leyend 
        xmin, xmax = ax1.get_xlim()
        ymin, ymax = ax1.get_ylim()
        x_pos = 0.6 * (xmin + xmax)
        y_pos = ymin + 0.69 * (ymax - ymin)
                
        # Add the metrics text (WSSR, DoF, SSR, R²) at dynamic position
        ax1.text(x_pos, y_pos,
                 f"R²={R2:.4f}",
                 fontsize=12)
        ##plot for multiple peaks during the fitting process (all the 17's peaks)
        
        grayscale = plt.cm.Greys(np.linspace(0.3, 0.3, n_peak))
        
        for iPeak in range(n_peak-1):
            if params[4 * iPeak + 1] > 1:  # should be bigger than 1, 
                peak = pearson7(wavelength, params[4 * iPeak:4 * iPeak + 4])
                #ax1.plot(wavelength, peak, ':', label=f'Peak {iPeak + 1}')
                ax1.plot(wavelength, peak, color=grayscale[iPeak], label=f'Peak {iPeak + 1}')
        ax1.set_title(f"{fl_name} {i_iter+1}/{n_iter} iteration", fontsize=13)
        
        # Residuals
        # === Bottom Subplot ===
        ax2 = plt.subplot(gs[1])
        ax2.plot(wavelength, intensity - curve, label='Residuals')
        ax2.axhline(0, color='black', linestyle='--')
        ax2.set_xlabel('Wavelength [nm]',fontsize=16)
        ax2.set_ylabel('Residuals [-]',fontsize=16)
        ax2.legend(loc=1)
        ax2.tick_params(labelsize=16)  # Larger font for x and y tick numbers
        ax2.get_legend().remove()
        plt.subplots_adjust(hspace=0.06)

        from scipy.integrate import trapz
        area_op = []
        for iPeak in range(n_peak-1): ###ipeak from zero
            if 8 <= iPeak <= 11:  # Calculate area for specific dye  peaks  (0-7)- NOMS. (8-11) DYES, (12-16) NOMS.
                peak = pearson7(wavelength, params[4 * iPeak:4 * iPeak + 4])
                area_op.append(trapz(peak, wavelength)) #### wavelength vs concentration
        if i_iter == n_iter - 1:
            print("Saving figure for the last iteration...")
            plt.savefig(f"{fig_dir}/{fl_name}_iteration_{i_iter + 1}.png", dpi=450)
        plt.show()
        plt.close()
    return area_op

############## END ######################################################