Includes code for simulation of photon detections, reflectivity estimation, and depth estimation. 

### Usage:  

*1. Load ground truth data as 2-D images*  
	- `I_truth`  
	- `D_truth` (in meters)  

*2. Generate simulated data*  
	See example code in `simulate_fdt.m` for setting simulation parameters.  

	`[k, t, detections, signal_mask] = simulate_detections(I_truth, D_truth, N, signal_level, noise_level, Tr, pulse_args);`  

	Inputs:  
	- `I_truth`: ground truth reflectivity image (2D matrix, `n1` x `n2`).  
	- `D_truth`: ground truth depth image in meters (`n1` x `n2`). If detection times are not needed (only doing reflectivity estimation), input depth as `[]`.  
	- `N`: number of laser pulses.  
	- `signal_level`: expected number of detected signal photons from unity-reflectivity pixel in one pulse repetition period, i.e. `q.*S`, where `q` is detection efficiency and `S` is average number of backreflected signal photons from unity-reflectivity pixel in one pulse repetition period. Should be either scalar or `n1` x `n2` matrix.  
	- `noise_level`: expected number of detected noise photons in one pulse repetition period (`q.*B + d` from paper). Should be either scalar or `n1` x `n2` matrix.  
	- `Tr`: pulse repetition period in seconds.  
	- `time_res`: time resolution of detector in seconds.  
	- `pulse_args`: two options...  
		1. provide pulse shape as histogram:  
			`pulse_args = {'pulse_hist', pulse_hist, pulse_bins};`  
			- `pulse_hist`: pulse shape as histogram.  
			- `pulse_bins`: time bins for pulse_hist.  
		2. provide parameters to simulate truncated Gaussian pulse shape:  
			`pulse_args = {'pulse_gauss', pulse_sigma, pulse_trunc};`  
			- `pulse_sigma`: RMS pulse width (standard deviation of Gaussian).  
			- `pulse_trunc`: minimum and maximum pulse times (bounds of truncated Gaussian), as `[pulse_min, pulse_max]`.  

	Outputs:  
	- `k`: number of detections at each pixel. (`n1` x `n2`)  
	- `t`: 3D array (`n1` x `n2` x `k_max`) with detection times as the first `k(i,j)` elements at each pixel `(i,j)`. The remaining elements are `NaN`.  
	- `detections`: sparse double (`n1*n2` x `N`) containing detection times.  
	- `signal_mask`: sparse logical (`n1*n2` x `N`), where ones indicate signal detections.  


*3. Estimate reflectivity*  
	Set parameters for SPIRAL algorithm:  
	- `tau`: strength of penalty/smoothing term  
	- `ainit`: initial value of alpha (see SPIRALTAP paper)  
	- `max_iter`: maximum number of iterations before terminating.  
	- `pen_type`: type of penalty/smoothing term `{'canonical','onb','rdp','rdp-ti','tv'}`. We use `'tv'`.  
	- `noisetype`: type of data likelihood term `{'poisson','binomial'}`. Reflectivity is estimated from number of detections at each pixel, which is a binomial distribution that can be approximated as a Poisson distribution.  

	Set observation matrix used for data likelihood term:  
	```
	A = @(x) signal_level.*x + noise_level; 
	AT = @(x) signal_level.*x; 
	```
	
	Call `SPIRALTAP_NEW` with parameters:  
	`I_MAP = SPIRALTAP_NEW(k,A,tau,N, ...optionalInputs...);`  

	See example code for usage of `SPIRALTAP_NEW`.  

*4. Censor noise detections using ROM filtering*

	Set a threshold value - this should be a (`n1` x `n2`) matrix that depends on the reflectivity at each pixel.

	Compute ROM image: `rom = get_rom(t, window_size);`

	Censor detections based on ROM and threshold: `t_signal = censor(t, rom, threshold);`

*5. Estimate depth*

	Use `t_signal_avg = nanmean(t_signal,3)` as input to SPIRALTAP, and use `rom` as initialization for SPIRALTAP. Because SPIRALTAP does not like `NaN` values, must first set all `NaN` values to `0`. 

	Set parameters for SPIRALTAP: (see descriptions for reflectivity estimation)  
	- `noisetype`: If pulse shape is approximately Gaussian, then noisetype should be `'gaussian'`.

	Simulation parameters do not contribute to time-of-arrival data likelihood term:  
	```
	A = @(x) x;
	AT = @(x) x;
	```

	Call `SPIRALTAP` with parameters:  
	`D_MAP = SPIRALTAP(t_signal_avg, A, tau, ...optionalInputs...);`  

	Convert from time bins to depth (meters):  
	`D_MAP = D_MAP * time_res * 3e8 * 0.5;`   