/**
 * bandlim.js
 * 
 * Compute bandwidth limited pulse duration from power spectrum.
 * Designed for visible regime, but should work elsewhere too.
 * 
 * @author: Marcelo alcocer
 * 
 */

// Module namepace
var bandlim = {

	// Speed of light [nm/fs]
	c : 299.792458,

	// Current Pulse instance
	pulse : null,
	
	// Charts
	chart_frequency: null,
	chart_time: null,

	// Generate FFT axes
	fftfreq : function(n, d){
		var f = null;
		if(n % 2 == 0){
			// Even length
			f = numeric.linspace(0, (n / 2) - 1);
			f = f.concat(numeric.linspace(-(n / 2), -1));
		}else{
			// Odd length
			f = numeric.linspace(0, (n - 1) / 2);
			f = f.concat(numeric.linspace(-(n - 1) / 2, -1));
		}
		f = numeric.div(f, (d * n));
		return(f);
	},

	// Shift FFT negative frequencies
	fftshift : function(x){

		function fftshift_(x_){
			var n = x_.length;

			// Even length
			if(n % 2 == 0){
				pos = x_.splice(0, (n / 2));
				x_ = x_.concat(pos);

				// Odd length
			}else{
				pos = x_.splice(0, ((n - 1) / 2) + 1);
				x_ = x_.concat(pos);
			}

			return(x_)
		}

		// numeric.T
		if(x instanceof numeric.T){
			x.x = fftshift_(x.x)
			x.y = fftshift_(x.y)

			// Array
		}else{
			x = fftshift_(x);
		}
		return(x);

	},

	// Pack arrays for plotting
	pack: function(x,y){
		var data = [];
		for(var i=0; i < x.length; i++){
			data.push([x[i],y[i]]);
		}
		return(data);
	},
	
	// Havercosine function
	hvc: function(x){
		return(0.5 * (1 + Math.cos(x)));
	},
	
	// Gaussian
	Gaussian : function(A, x0, sigma){

		// Amplitude
		this.A = null;
		
		// Centre
		this.x0 = null
		
		// Standard deviation
		this.sigma = null;

		// Dependent variable
		this.y = function(x){
			var y = [];
			for(var i = 0; i < x.length; i++){
				y.push(
					this.A *
					Math.exp(
						-Math.pow(x[i] - this.x0, 2) / (2 * Math.pow(this.sigma, 2))
					)
				);
			}
			return(y);
		};
		
		// Standard deviation function
		this.fwhm2sigma = function(fwhm_){
			return(fwhm_ / (2 * Math.sqrt(2 * Math.LN2)));
		};
		
		// Full width at half maximum functino
		this.sigma2fwhm = function(sigma_){
			return(2 * Math.sqrt(2 * Math.LN2) * sigma_);
		};

		// Full width at half maximum
		Object.defineProperties(this, {
			"fwhm": {
				"get": function(){
					return(this.sigma2fwhm(this.sigma))
				}
			}
		});
		
		// Init
		this.init = function(){
			
			// Assign properties from argument
			this.A = A;
			this.x0 = x0;
			this.sigma = sigma;
			
		};
		this.init();

	},
	
	// Tukey window
	Tukey: function(alpha){
		
		// Alpha
		this.alpha = null;
		
		// Dependent variable
		this.y = function(N){
			var y = [];
			for(var i = 0; i < N; i++){
				if(i < ((this.alpha * (N-1))/2)){
					y.push(
						bandlim.hvc(
							Math.PI * (
								((2*i)/(this.alpha*(N-1))) - 1
							)
						)
					);
				}else if(i <= ((N-1)*(1-(this.alpha/2)))){
					y.push(1);
				}else{
					y.push(
						bandlim.hvc(
							Math.PI * (
								((2*i)/(this.alpha*(N-1))) - (2/this.alpha) + 1
							)
						)
					);
				}
			}
			return(y);
		};
		
		// Init
		this.init = function(){

			// Assign properties from argument
			this.alpha = alpha;
			
		};
		this.init();
		
	},
	
	// General spectral/tmeporal profile
	Profile: function(x, y){

		// Independent axis (array)
		this.axis = null;
		
		// Intensity (array)
		this.intensity = null;
		
		// Full-width at half-maximum
		Object.defineProperties(this, {
			"fwhm" : {
				"get" : function(){

					// Locate FWHM indices
					var y = this.intensity;
					var max = Math.max.apply(null, y);
					var inds = {
						max: y.indexOf(max),
						low: null,
						high: null
					};
					inds.low = y.slice(0, inds.max).findIndex(function(val){
						return(val > (max / 2))
					});
					inds.high = inds.max + y.slice(inds.max).findIndex(function(val){
						return(val < (max / 2))
					});
					
					// Interpolate
					var margin = 2, low, high;
					try{
						low = numeric.spline(
							y.slice(inds.low - margin, inds.low + margin),
							this.axis.slice(inds.low - margin, inds.low + margin)
						).at(max/2);
						high = numeric.spline(
							y.slice(inds.high - margin, inds.high + margin),
							this.axis.slice(inds.high - margin, inds.high + margin)
						).at(max/2);
					}catch(e){
						low = Number.NaN;
						high = Number.NaN;
					}
					
					// Return
					return([low, high]);
				}
			}
		});
		
		// Fit
		this.fit = function(f, x0){

			// Define sum squared error objective function 
			var this_ = this;
			obj = function(x_){
				return(
					numeric.sub(this_.intensity, f(this_.axis, x_)).reduce(
						function(total, val){
							return(total + Math.pow(val,2));
						}
					)
				);
			}
			
			// (Unconstrained) optimisation
			var tol = 1e-8;
			var maxiter = 1000;
			var x = numeric.uncmin(obj, x0, tol, undefined, maxiter);
			
			// Return fit object
			return(x);
			
		},

		// Plot
		this.plot = function(chart, id, type="line"){
			
			// Construct series options
			var opt_series = {
				id: id,
				name: id,
				className: id.replace(" ",""),
				type: type,
				animation: false,
				data: bandlim.pack(
					this.axis,
					this.intensity
				)
			};
			
			// Create new series
			var series = chart.get(id);
			if(typeof series == "undefined"){
				chart.addSeries(opt_series);
				
			// Update existing series
			}else{
				series.update(opt_series)
			}
		};

		
		// Init
		this.init = function(){
			
			// Assign properties from constructor arguments
			this.axis = x;
			this.intensity = y;

		};
		this.init();
		
	},

	// Window
	Window: function(func, x, limits){
		
		// Window function
		this.function_ = null;
		
		// Window limits
		this._limits = null;
		Object.defineProperties(this,{
			"limits": {
				"get": function(){return(this._limits)},
				"set": function(limits_){
					
					// Sort and set private limits
					limits_.sort(function(a, b){
						return(a - b)
					});
					this._limits = limits_;

					// Linearise axis
					var x = this.axis;
					var x_lin = numeric.linspace(
						x[0],
						x[x.length - 1],
						x.length
					);
					
					// Find window limits
					var inds = [
						x_lin.findIndex(function(val){
							return(val >= limits_[0]);
						}),
						x_lin.findIndex(function(val){
							return(val >= limits_[1]);
						})
					];
					
					// Generate window within limits
					var y_lin = [];
					var i = x.length;
					while(i--) y_lin[i] = 0;
					var win = this.function_(inds[1] - inds[0] + 1);
					for(i = 0; i < win.length; i++){
						y_lin[inds[0] + i] = win[i];
					}
					
					// Set window values
					this.intensity = numeric.spline(
						x_lin,
						y_lin
					).at(x)

				}
			}
		});

		// Init
		this.init = function(){
			
			// Assign properties from arguments
			this.function_ = func;
			this.axis = x;			// Needed for comuptation in limits setter
			this.limits = limits;
			
			// Call superclass constructor with interpolated window
			bandlim.Profile.call(
				this,
				this.axis,
				this.intensity
			);
			
		};
		this.init();
		
	},
	
	// Electric Field
	Field: function(x, y, from_amplitude=false){

		// Amplitude
		this._amplitude = null;
		Object.defineProperties(this,{
			"amplitude": {
				"get": function(){return(this._amplitude)},
				"set": function(amplitude_){

					// Set private property
					this._amplitude = amplitude_;
					
					// Compute intensity
					this._intensity = this.amplitude.mul(this.amplitude.conj()).x;

				}
			}
		});
		
		// Intensity
		this._intensity = null;
		Object.defineProperties(this,{
			"intensity": {
				"get": function(){return(this._intensity)},
				"set": function(intensity_){
					if(intensity_ != null){
					
						// Truncate and set private property
						this._intensity = intensity_.map(
							function(val){
								return(val < 0 ? 0 : val);
							}
						);
						
						// Compute amplitude
						this._amplitude = new numeric.T(
							numeric.sqrt(this.intensity),
							numeric.rep([this.intensity.length], 0)
						);

					}
				}
			}
		});
		
		// Init
		this.init = function(){
			
			// Assign properties from arguments
			if(from_amplitude){
				this.amplitude = y;
			}else{
				this.intensity = y;
			}
			
			// Assign intensity (superclass constructor)
			bandlim.Profile.call(this, x, this.intensity);
			
		};
		this.init();
		
	},
	
	// Ultrafast pulse
	Pulse: function(frequency, intensity){
		
		// Spectrum
		this.spectrum = null;
		
		// Spectral window
		this.window = null;
		
		// Electric field — frequency domain
		this.field_frequency = null;
		
		// Electric field — time domain
		Object.defineProperties(this, {
			"field_time" : {
				"get" : function(){

					// Zero frequency axis at start of window
					var win = this.window;
					var dlim = win.limits[1] - win.limits[0];
					var f = numeric.sub(
						this.field_frequency.axis,
						win.limits[0]
					);

					// Define linear conjugate axes
					var t_ny = 100, dt = 0.1
					var t_lin = numeric.linspace(
						-t_ny,
						t_ny,
						((2 * t_ny) / dt) + 1
					);
					var f_lin = bandlim.fftshift(
						bandlim.fftfreq(
							t_lin.length,
							t_lin[1] - t_lin[0]
						)
					);
					numeric.subeq(f_lin, f_lin[0]) 		// Start at 0

					// Interpolate to linear frequency axis
					var ind = f_lin.findIndex(function(val){
						return(val > dlim);
					});
					var A = this.field_frequency.amplitude;
					var A_lin = {
						x : [],
						y : []
					};
					A_lin.x = numeric.spline(f, A.x).at(f_lin.slice(0, ind));
					A_lin.y = numeric.spline(f, A.y).at(f_lin.slice(0, ind));
					
					// Zero pad
					for(var i = ind; i < f_lin.length; i++){
						A_lin.x.push(0);
						A_lin.y.push(0);
					}
					A_lin = new numeric.T(A_lin.x, A_lin.y);

					// Transform
					return(
						new bandlim.Field(
							t_lin,
							bandlim.fftshift(A_lin.fft()),
							true
						)
					);

				}
			}
		});

		// Fits
		this.fit_frequency = null;
		this.fit_time = null;
		
		// Redraw
		this.redraw= function(){
			
			// Transform
			var field_time = this.field_time;
			
			// Update FWHMs
			var fwhm = {
				limits: {
					frequency: this.field_frequency.fwhm,
					time: field_time.fwhm
				},
				value: {
					frequency: null,
					time: null
				}
			};
			fwhm.value.frequency = fwhm.limits.frequency[1] - fwhm.limits.frequency[0];
			fwhm.value.time= fwhm.limits.time[1] - fwhm.limits.time[0];
			
			// Update fits
			var fit = {
				frequency: new bandlim.Gaussian(),
				time: new bandlim.Gaussian()
			};
			this.field_frequency.fit(
				function(x, params){
					fit.frequency.A = params[0];
					fit.frequency.x0 = params[1];
					fit.frequency.sigma = params[2];
					return(fit.frequency.y(x));
				},
				[
					1.0,
					(this.window.limits[0] + this.window.limits[1])/2,
					fit.frequency.fwhm2sigma(fwhm.value.frequency)
				]
			);
			field_time.fit(
				function(x, params){
					fit.time.A = params[0];
					fit.time.x0 = params[1];
					fit.time.sigma = params[2];
					return(fit.time.y(x));
				},
				[
					Math.max.apply(null, field_time.intensity),
					0,
					fit.time.fwhm2sigma(fwhm.value.time)
				]
			);
			
			// Update FWHM and TBWP outputs
			var digits = 3;
			document.getElementById("fwhm_frequency_numerical").innerHTML = fwhm.value.frequency.toFixed(digits);
			document.getElementById("fwhm_time_numerical").innerHTML = fwhm.value.time.toFixed(digits);
			document.getElementById("tbwp_numerical").innerHTML = (fwhm.value.time * fwhm.value.frequency).toFixed(digits);
			document.getElementById("fwhm_frequency_fit").innerHTML = fit.frequency.fwhm.toFixed(digits);
			document.getElementById("fwhm_time_fit").innerHTML = fit.time.fwhm.toFixed(digits);
			document.getElementById("tbwp_fit").innerHTML = (fit.time.fwhm * fit.frequency.fwhm).toFixed(digits);
			
			// Update series
			this.window.plot( bandlim.chart_frequency, "Window");
			this.field_frequency.plot(
				bandlim.chart_frequency,
				"Field Intensity",
				"area"
			);
			field_time.plot(
				bandlim.chart_time,
				"Field Intensity",
				"area"
			);
			(new bandlim.Profile(
				this.field_frequency.axis,
				fit.frequency.y(this.field_frequency.axis)
			)).plot(bandlim.chart_frequency, "Fit");
			(new bandlim.Profile(
				field_time.axis,
				fit.time.y(field_time.axis)
			)).plot(bandlim.chart_time, "Fit");
			
			// Update plot bands
			var id = {
				plotBand: {
					numerical: "PlotBand Numerical",
					fit: "PlotBand Fit",
				},
				label: {
					numerical: "Label Numerical",
					fit: "Label Fit",
				}
			};
			bandlim.chart_frequency.xAxis[0].update(
				{
					plotBands: [{
						id: id.plotBand.numerical,
						className: id.plotBand.numerical.replace(" ",""),
						from: fwhm.limits.frequency[0],
						to: fwhm.limits.frequency[1],
						label: {
							align: "left",
							textAlign: "right",
							x: -5,
							id: id.label.numerical,
							className: id.label.numerical.replace(" ",""),
							text: "Numerical<br/>" + fwhm.value.frequency.toFixed(digits) + " PHz"
						}
					},{
						id: id.plotBand.fit,
						className: id.plotBand.fit.replace(" ",""),
						from: fit.frequency.x0 - (fit.frequency.fwhm/2),
						to: fit.frequency.x0 + (fit.frequency.fwhm/2),
						label: {
							align: "right",
							textAlign: "left",
							x: 5,
							id: id.label.fit,
							className: id.label.fit.replace(" ",""),
							text: "Fit<br/>" + fit.frequency.fwhm.toFixed(digits) + " PHz"
						}
					}]
				}
			);
			bandlim.chart_time.xAxis[0].update(
				{
					plotBands: [{
						id: id.plotBand.numerical,
						className: id.plotBand.numerical.replace(" ",""),
						from: fwhm.limits.time[0],
						to: fwhm.limits.time[1],
						label: {
							align: "left",
							textAlign: "right",
							x: -5,
							id: id.label.numerical,
							className: id.label.numerical.replace(" ",""),
							text: "Numerical<br/>" + fwhm.value.time.toFixed(digits) + " fs"
						}
					},
					{
						id: id.plotBand.fit,
						className: id.plotBand.fit.replace(" ",""),
						from: fit.time.x0 - (fit.time.fwhm/2),
						to: fit.time.x0 + (fit.time.fwhm/2),
						label: {
							align: "right",
							textAlign: "left",
							x: 5,
							id: id.label.fit,
							className: id.label.fit.replace(" ",""),
							text: "Fit<br/>" + fit.time.fwhm.toFixed(digits) + " fs"
						}
					}]
				}
			);
			
		};
		
		// Init
		this.init = function(){
			
			// Instantiate spectrum
			this.spectrum = new bandlim.Profile(
				frequency,
				numeric.div(
					intensity,
					Math.max.apply(null, intensity)
				)
			);
			this.spectrum.plot(
				bandlim.chart_frequency,
				"Spectrum"
			);
			
			// Instantiate window
			var tukey = new bandlim.Tukey(0.5);
			this.window = new bandlim.Window(
				function(N){
					return(tukey.y(N));
				},
				this.spectrum.axis,
				[
					frequency[0],
					frequency[frequency.length - 1]
				]
			);
			
			// Instantiate (windowed) electric field — frequency domain
			this.field_frequency = new bandlim.Field(
				this.spectrum.axis,
				numeric.mul(
					this.spectrum.intensity,
					this.window.intensity
				),
				false
			);
			
			// Redraw
			this.redraw();
			
		};
		this.init();
		
	},

	// AvaSoft file
	AvaSoftFile : function(buffer){

		this.start_pixel = null;

		this.stop_pixel = null;

		this.wavelength = null;

		this.sample = null;

		this.dark = null;

		this.reference = null;

		Object.defineProperties(this, {
			"frequency" : {
				"get" : function(){
					return(numeric.div(bandlim.c, this.wavelength))
				}
			},
			"irradiance" : {
				"get" : function(){
					return(numeric.div(numeric.sub(this.sample, this.dark),
							this.reference))
				}
			}
		});

		// Initialise from ArrayBuffer
		this.init = function(){
			this.magic = new Uint8Array(buffer.slice(0,3));
			this.start_pixel = (new Uint16Array(buffer.slice(89, 91)))[0];
			this.stop_pixel = (new Uint16Array(buffer.slice(91, 93)))[0];
			var sz_float32 = 4;
			var pixel_bytes = (this.stop_pixel - this.start_pixel + 1)
					* sz_float32;
			var data_start = 328;
			this.wavelength = new Float32Array(buffer.slice(data_start,
					data_start + pixel_bytes));
			this.sample = new Float32Array(buffer.slice(data_start
					+ pixel_bytes, data_start + (2 * pixel_bytes)));
			this.dark = new Float32Array(buffer.slice(data_start
					+ (2 * pixel_bytes), data_start + (3 * pixel_bytes)));
			this.reference = new Float32Array(buffer.slice(data_start
					+ (3 * pixel_bytes), data_start + (4 * pixel_bytes)));
		};
		this.init();

	},

	// Module init
	init : function(){

		// Create charts
		var opt_tooltip = {
			enabled : false
		};
		var opt_axis = {
			type : "linear"
		} ;
		bandlim.chart_frequency = Highcharts.chart(
			"chart_frequency",
			{
				chart : {
					type : "line",
					zoomType : "xy"
				}, 
				tooltip : opt_tooltip,
				xAxis : [
					{
						type: "linear",
						title : {text : "Frequency [PHz]"}
					},
					{
						type: "linear",
						title : {text : "Wavelength [nm]"},
						opposite: true,
						linkedTo: 0,
						tickPositioner: function(){
							var positions = [200,300,400,500,600,700,800,900,1000,1100];
							return(positions.map(function(val){
								return(bandlim.c / val);
							}));
						},
						labels:{
							formatter: function(){
								return(
									(bandlim.c / Number(this.value)).toFixed(0)
								);
							}
						}
					}
				],
				yAxis : {
					type: "linear",
					title : {text : "Intensity [a.u.]"}
				},
				title: {
					text : "Frequency Domain"
				}
			}
		);
		bandlim.chart_time = Highcharts.chart(
			"chart_time",
			{
				chart : {
					type : "line",
					zoomType : "xy"
				}, 
				tooltip : opt_tooltip,
				xAxis : {
					type: "linear",
					title : {text : "Time [fs]"}
				},
				yAxis : {
					type: "linear",
					title : {text : "Intensity [a.u.]"}
				},
				title: {
					text : "Time Domain"
				}
			}
		);
		
		// File handling
		document.getElementById("file").addEventListener(
			"change",
			function(changeEvent){

				
				// Determine file type
				var file = changeEvent.target.files[0];
				var reader = new FileReader();
				var onloadend = function(loadEvent){};
				var read = function(){};
				if(file != undefined){
					switch(file.type){

						// CSV
						case "text/csv":
						case "text/plain":
							
							// Set onloadend callback
							reader.onloadend = function(loadEvent){
								if(loadEvent.target.readyState == FileReader.DONE){
									
									// Extract spectrum
									var csv = numeric.parseCSV(loadEvent.target.result);
									csv.sort(function(a, b){
										if(a[0] === b[0]){
											return(0);
										}else{
											return(a[0] > b[0] ? -1 : 1);
										}
									});
									var frequency = [], intensity = [];
									for(var i=0; i < csv.length; i++){
										frequency.push(
											bandlim.c / csv[i][0]
										);
										intensity.push(csv[i][1]);
									}
									
									// Create pulse instance
									bandlim.pulse = new bandlim.Pulse(
										frequency,
										intensity
									)

								}
							};
							
							// Read file
							read = reader.readAsText(file);
							break;

						// AvaSoft 8
						case "":
						case "application/octet-stream":
							
							// Set onloadend callback
							reader.onloadend = function(loadEvent){
								if(loadEvent.target.readyState == FileReader.DONE){

									// Extract spectrum
									var avs_spec = new bandlim.AvaSoftFile(
										loadEvent.target.result
									);
									
									// Create Pulse instance
									bandlim.pulse = new bandlim.Pulse(
										avs_spec.frequency.reverse(),
										avs_spec.irradiance.reverse()
									);
								}
							};

							// Read file
							reader.readAsArrayBuffer(file);
							break;

						// Unsupported
						default:
							throw("Invalid file type");
					}
				}

			},
			false
		);
		

		// Window limits
		function updateWindowFactory(j){
			return(
				function(changeEvent){
					
					// Update window limits
					bandlim.pulse.window.limits = [
						Number(changeEvent.target.value),
						bandlim.pulse.window.limits[(j + 1) % 2]
					];
					
					// Update frequency domain electric field
					bandlim.pulse.field_frequency.intensity = numeric.mul(
						bandlim.pulse.spectrum.intensity,
						bandlim.pulse.window.intensity
					);
					
					// Redraw
					bandlim.pulse.redraw();

				}
			)
		}
		for(var i=0; i < 2; i++){
			document.getElementById("limits_" + i).addEventListener(
				"change",
				updateWindowFactory(i),
				false
			);
		}
		
		// Initial demo spectrum
		var f = [0.6200079623, 0.6195851073, 0.6191500413, 0.6187283553, 0.6183072433, 0.6178867042, 0.6174667367, 0.6170346396, 0.6166158292, 0.6161975869, 0.6157799117, 0.6153628022, 0.6149336437, 0.6145176796, 0.6141022778, 0.6136874373, 0.6132731569, 0.6128469071, 0.6124337603, 0.6120211702, 0.6116091356, 0.6111976555, 0.6107867286, 0.6103639269, 0.60995412, 0.6095448631, 0.6091361549, 0.6087279945, 0.6083203807, 0.6079009855, 0.6074944781, 0.607088514, 0.6066830922, 0.6062782115, 0.6058738708, 0.6054578411, 0.6050545935, 0.6046518828, 0.6042497077, 0.6038480674, 0.6034469605, 0.6030342559, 0.6026342292, 0.6022347328, 0.6018357658, 0.601437327, 0.6010394154, 0.60064203, 0.6002331518, 0.5998368315, 0.5994410342, 0.5990457588, 0.5986510044, 0.59825677, 0.5978630544, 0.5974579497, 0.5970652845, 0.5966731351, 0.5962815005, 0.5958903796, 0.5954997716, 0.5951096752, 0.594708292, 0.5943192316, 0.5939306799, 0.593542636, 0.5931550987, 0.5927680672, 0.5923815405, 0.5919955175, 0.5916099972, 0.5912133194, 0.5908288169, 0.5904448143, 0.5900613104, 0.5896783045, 0.5892957954, 0.5889137823, 0.5885322641, 0.5881512399, 0.5877591852, 0.587379161, 0.586999628, 0.5866205851, 0.5862420314, 0.5858639659, 0.5854863878, 0.585109296, 0.5847326897, 0.5843565678, 0.5839695541, 0.5835944131, 0.5832197537, 0.5828455751, 0.5824718763, 0.5820986564, 0.5817259144, 0.5813536496, 0.5809818608, 0.5806105473, 0.5802397081, 0.5798581265, 0.5794882476, 0.5791188402, 0.5787499035, 0.5783814365, 0.5780134385, 0.5776459084, 0.5772788454, 0.5769122486, 0.5765461172, 0.5761804501, 0.5758152466, 0.5754505058, 0.5750862267, 0.574711391, 0.5743480468, 0.5739851618, 0.573622735, 0.5732607656, 0.5728992528, 0.5725381956, 0.5721775933, 0.5718174449, 0.5714577496, 0.5710985065, 0.5707397148, 0.5703813737, 0.5700234822, 0.5696660396, 0.569309045, 0.5689524975, 0.5685963964, 0.5682407407, 0.5678855297, 0.5675307625, 0.5671657081, 0.5668118394, 0.5664584121, 0.5661054252, 0.5657528779, 0.5654007695, 0.5650490991, 0.5646978658, 0.564347069, 0.5639967077, 0.5636467812, 0.5632972887, 0.5629482292, 0.5625996022, 0.5622514066, 0.5619036418, 0.5615563069, 0.5612094011, 0.5608629237, 0.5605168739, 0.5601712508, 0.5598260537, 0.5594812817, 0.5591369342, 0.5587930103, 0.5584495092, 0.5581064301, 0.5577637723, 0.5574215351, 0.5570797176, 0.556738319, 0.5563973386, 0.5560567756, 0.5557166293, 0.5553768989, 0.5550375835, 0.5546986826, 0.5543601953, 0.5540221208, 0.5536844584, 0.5533472074, 0.5530103669, 0.5526739363, 0.5523379148, 0.5520023016, 0.551667096, 0.5513322973, 0.5509979048, 0.5506639176, 0.550330335, 0.5499971564, 0.5496643809, 0.5493320079, 0.5490000366, 0.5486684663, 0.5483372963, 0.5480065258, 0.5476861595, 0.5473561739, 0.5470265857, 0.5466973941, 0.5463685985, 0.5460401982, 0.5457121924, 0.5453845804, 0.5450573616, 0.5447305351, 0.5444041004, 0.5440780567, 0.5437524032, 0.5434271394, 0.5431022645, 0.5427777778, 0.5424536786, 0.5421299662, 0.5418066399, 0.5414836991, 0.541161143, 0.5408487281, 0.5405269279, 0.5402055103, 0.5398844748, 0.5395638206, 0.5392435471, 0.5389236536, 0.5386041393, 0.5382850038, 0.5379662462, 0.5376478659, 0.5373298622, 0.5370122344, 0.536694982, 0.5363877011, 0.5360711859, 0.5357550441, 0.5354392749, 0.5351238777, 0.5348088519, 0.5344941967, 0.5341799116, 0.5338659959, 0.5335524489, 0.53323927, 0.5329359323, 0.5326234765, 0.5323113869, 0.5319996628, 0.5316883036, 0.5313773087, 0.5310666773, 0.5307564089, 0.5304465028, 0.5301369584, 0.529837139, 0.5295283052, 0.5292198312, 0.5289117164, 0.5286039602, 0.5282965619, 0.527989521, 0.5276828367, 0.5273765085, 0.5270798024, 0.5267741737, 0.5264688993, 0.5261639784, 0.5258594106, 0.5255551952, 0.5252513316, 0.5249478191, 0.5246538387, 0.5243510162, 0.524048543, 0.5237464186, 0.5234446423, 0.5231432136, 0.5228421319, 0.5225505046, 0.5222501045, 0.5219500496, 0.5216503393, 0.521350973, 0.5210519501, 0.52075327, 0.5204639676, 0.5201659611, 0.5198682956, 0.5195709705, 0.5192739854, 0.5189773396, 0.5186810325, 0.5183940274, 0.5180983859, 0.5178030813, 0.5175081132, 0.517213481, 0.5169191841, 0.5166341249, 0.5163404867, 0.5160471822, 0.5157542106, 0.5154615715, 0.5151692643, 0.5148861314, 0.5145944763, 0.5143031514, 0.5140121562, 0.5137214901, 0.5134311526, 0.5131499264, 0.5128602344, 0.5125708692, 0.5122818305, 0.5119931175, 0.511713464, 0.5114253911, 0.5111376424, 0.5108502173, 0.5105631152, 0.5102763357, 0.509998554, 0.5097124082, 0.5094265833, 0.5091410788, 0.5088558941, 0.5085796563, 0.5082951, 0.508010862, 0.5077269417, 0.5074433386, 0.5071686319, 0.506885652, 0.5066029876, 0.5063206384, 0.5060471456, 0.5057654154, 0.5054839988, 0.5052028951, 0.504922104, 0.5046501195, 0.5043699423, 0.504090076, 0.5038105201, 0.5035312741, 0.5032607856, 0.5029821485, 0.5027038198, 0.502425799, 0.5021564965, 0.5018790806, 0.501601971, 0.5013251672, 0.5010570431, 0.5007808402, 0.5005049417, 0.5002293471, 0.4999623935, 0.4996873958, 0.4994127005, 0.4991383071, 0.498864215, 0.4985987161, 0.4983252161, 0.498052016, 0.4977791153, 0.4975147698, 0.4972424574, 0.4969704429, 0.4967069554, 0.4964355263, 0.4961643938, 0.4958935572, 0.49563121, 0.4953609551, 0.4950909947, 0.4948213284, 0.4945601142, 0.4942910257, 0.4940222299, 0.4937537263, 0.4934936378, 0.4932257083, 0.4929580696, 0.4926988183, 0.492431751, 0.492164973, 0.4918984839, 0.4916403457, 0.4913744243, 0.4911087904, 0.4908434435, 0.4905864112, 0.4903216283, 0.4900571312, 0.4898009215, 0.4895369856, 0.4892733341, 0.4890179431, 0.4887548502, 0.4884920402, 0.4882295127, 0.48797521, 0.4877132376, 0.4874515463, 0.4871980531, 0.4869369143, 0.4866760552, 0.4864233677, 0.4861630585, 0.4859030277, 0.4856432749, 0.4853916584, 0.4851324519, 0.4848735221, 0.4846227025, 0.4843643164, 0.4841062057, 0.4838561791, 0.4835986095, 0.483341314, 0.4830920766, 0.4828353197, 0.4825788355, 0.4823303837, 0.4820744356, 0.4818187589, 0.4815710889, 0.4813159458, 0.4810610729, 0.480814181, 0.4805598391, 0.4803057661, 0.4800596487, 0.4798061042, 0.4795528273, 0.4793074808, 0.4790547299, 0.4788022455, 0.4785576662, 0.4783057054, 0.4780540097, 0.477810194, 0.4775590194, 0.4773157082, 0.4770650531, 0.4768146611, 0.4765721076, 0.4763222326, 0.4760726196, 0.4758308203, 0.4755817218, 0.4753328841, 0.4750918354, 0.4748435099, 0.4746029572, 0.4743551424, 0.4741075862, 0.4738677784, 0.4736207305, 0.4733739401, 0.4731348737, 0.4728885892, 0.4726500126, 0.4724042325, 0.4721587079, 0.471920867, 0.4716758445, 0.4714310762, 0.4711939677, 0.4709496992, 0.4707130745, 0.4704693042, 0.4702257862, 0.4699898882, 0.4697468662, 0.4695114483, 0.4692689207, 0.4690266435, 0.4687919468, 0.4685501618, 0.4683159416, 0.4680746471, 0.4678336012, 0.4676000967, 0.4673595392, 0.4671265075, 0.4668864369, 0.4666538767, 0.4664142915, 0.4661749522, 0.4659431, 0.4657042439, 0.4654728597, 0.4652344853, 0.464996355, 0.4647656734, 0.4645280227, 0.4642978054, 0.4640606328, 0.4638308785, 0.4635941825, 0.463357728, 0.463128669, 0.462892689, 0.4626640894, 0.4624285824, 0.4622004409, 0.4619654057, 0.4617377208, 0.4615031558, 0.461268829, 0.4610418301, 0.4608079713, 0.4605814257, 0.4603480337, 0.46012194, 0.4598890133, 0.4596633701, 0.4594309074, 0.4592057134, 0.4589737132, 0.4587419473, 0.4585174281, 0.4582861227, 0.4580620493, 0.4578312029, 0.4576075741, 0.4573771855, 0.4571539998, 0.4569240676, 0.4567013238, 0.4564718466, 0.4562495434, 0.4560205199, 0.455798656, 0.4555700848, 0.4553486588, 0.4551205386, 0.4548926469, 0.4546718788, 0.454444436, 0.4542241027, 0.4539971076, 0.4537772076, 0.4535506589, 0.4533311911, 0.4531050874, 0.4528860505, 0.4526603905, 0.4524417833, 0.4522165656, 0.4519983868, 0.4517736102, 0.4515558585, 0.4513315217, 0.4511141959, 0.4508902976, 0.4506733964, 0.4504499354, 0.4502334575, 0.4500104325, 0.4497943767, 0.4495785283, 0.4493561514, 0.4491407233, 0.4489187793, 0.4487037702, 0.4484822577, 0.4482676665, 0.4480465843, 0.4478324097, 0.4476117565, 0.4473979973, 0.4471777718, 0.4469644268, 0.4467446279, 0.4465316959, 0.4463123223, 0.4460998021, 0.4458874842, 0.4456687429, 0.4454568351, 0.445238516, 0.445027017, 0.444809119, 0.4445980276, 0.4443805493, 0.4441698644, 0.4439528048, 0.4437425252, 0.4435324447, 0.4433160074, 0.4431063305, 0.4428903088, 0.4426810342, 0.4424654269, 0.4422565536, 0.4420413595, 0.4418328863, 0.4416246096, 0.4414100299, 0.4412021516, 0.4409879821, 0.4407805011, 0.4405667406, 0.4403596557, 0.4401527653, 0.4399396132, 0.4397331172, 0.4395203712, 0.4393142685, 0.4391083591, 0.438896217, 0.4386906992, 0.4384789604, 0.438273833, 0.4380688975, 0.4378577584, 0.4376532117, 0.437442473, 0.437238314, 0.4370279746, 0.4368242022, 0.4366206198, 0.4364108742, 0.4362076768, 0.4359983275, 0.435795514, 0.4355928891, 0.4353841294, 0.4351818868, 0.434979832, 0.4347716594, 0.4345699852, 0.4343622046, 0.4341609101, 0.433959802, 0.4337526043, 0.4335518742, 0.4333450658, 0.4331447127, 0.4329445447, 0.4327383152, 0.4325385226, 0.4323389144, 0.4321332613, 0.4319340268, 0.431734976, 0.4315298969, 0.4313312183, 0.4311265226, 0.4309282152, 0.4307300901, 0.4305259643, 0.4303282089, 0.430130635, 0.4299270769, 0.4297298711, 0.4295328462, 0.4293298534, 0.429133195, 0.4289367166, 0.4287342867, 0.4285381734, 0.4283422395, 0.4281403702, 0.4279447997, 0.4277494079, 0.4275480968, 0.427353067, 0.427158215, 0.4269574598, 0.4267629683, 0.426568654, 0.4263684525, 0.4261744971, 0.4259807181, 0.425781068, 0.4255876466, 0.4253944008, 0.4251952998, 0.42500241, 0.4248096952, 0.424611141, 0.4244187808, 0.4242265948, 0.4240285852, 0.4238367523, 0.4236450929, 0.4234476256, 0.4232563179, 0.423065183, 0.4228742207, 0.4226774712, 0.4224868586, 0.4222964179, 0.4221002056, 0.4219101132, 0.4217201919, 0.4215245146, 0.4213349402, 0.4211455363, 0.4209563026, 0.4207613333, 0.4205724447, 0.4203837255, 0.4201892862, 0.4200009106, 0.4198127039, 0.4196246658, 0.4194309279, 0.4192432315, 0.4190557031, 0.4188683424, 0.418675302, 0.4184882812, 0.4183014274, 0.4181089091, 0.4179223938, 0.4177360449, 0.4175498621, 0.4173580348, 0.4171721887, 0.4169865081, 0.4168009927, 0.4166098527, 0.4164246722, 0.4162396562, 0.4160548046, 0.4158643482, 0.4156798297, 0.4154954749, 0.4153112835, 0.4151215071, 0.4149376471, 0.4147539498, 0.4145704151, 0.4143813151, 0.4141981099, 0.4140150668, 0.4138321853, 0.4136437579, 0.4134612043, 0.4132788117, 0.4130965799, 0.4129145088, 0.412726916, 0.4125451706, 0.4123635851, 0.4121821594, 0.4119952313, 0.4118141295, 0.4116331869, 0.4114524032, 0.4112717782, 0.4110856747, 0.4109053715, 0.4107252264, 0.4105452392, 0.4103654096, 0.4101801253, 0.4100006154, 0.4098212626, 0.4096420666, 0.4094630272, 0.4092785567, 0.4090998349, 0.4089212691, 0.4087428591, 0.4085646047, 0.4083809427, 0.4082030037, 0.4080252198, 0.4078475907, 0.4076701161, 0.4074872572, 0.4073100961, 0.4071330889, 0.4069562355, 0.4067795357, 0.4065974746, 0.4064210862, 0.4062448507, 0.4060687679, 0.4058928378, 0.40571706, 0.4055359486, 0.4053604797, 0.4051851627, 0.4050099972, 0.4048349831, 0.4046601201, 0.4044799509, 0.4043053945, 0.4041309886, 0.4039567332, 0.403782628, 0.4036086728, 0.4034348674, 0.4032557874, 0.4030822857, 0.4029089333, 0.4027357299, 0.4025626754, 0.4023897695, 0.4022170121, 0.4020390114, 0.401866555, 0.4016942464, 0.4015220856, 0.4013500723, 0.4011782063, 0.4010064874, 0.4008295562, 0.4006581356, 0.4004868616, 0.400315734, 0.4001447525, 0.399973917, 0.3998032273, 0.3996326833, 0.399456962, 0.3992867132, 0.3991166094, 0.3989466505, 0.3987768363, 0.3986071666, 0.3984376412, 0.39826826, 0.3980990227, 0.3979246473, 0.3977557018, 0.3975868997, 0.3974182409, 0.397249725, 0.3970813521, 0.3969131218, 0.3967450339, 0.3965770884, 0.396409285, 0.3962363865, 0.3960688712, 0.3959014975, 0.3957342653, 0.3955671742, 0.3954002242, 0.3952334151, 0.3950667466, 0.3949002187, 0.394733831, 0.3945675836, 0.3944014761, 0.3942355084, 0.3940645004, 0.3938988162, 0.3937332712, 0.3935678653, 0.3934025983, 0.39323747, 0.3930724803, 0.392907629, 0.3927429159, 0.3925783409, 0.3924139037, 0.3922496042, 0.3920854423, 0.3919214176, 0.3917575302, 0.3915937798, 0.3914250555, 0.3912615828, 0.3910982467, 0.3909350468, 0.3907719831, 0.3906090554, 0.3904462634, 0.3902836072, 0.3901210863, 0.3899587008, 0.3897964504, 0.389634335, 0.3894723543, 0.3893105083, 0.3891487967, 0.3889872194, 0.3888257762, 0.388664467, 0.3885032916, 0.3883422498, 0.3881813414, 0.3880205664, 0.3878599244, 0.3876994155, 0.3875390393, 0.3873787957, 0.3872186846, 0.3870587058, 0.3868988591, 0.3867391445, 0.3865795616, 0.3864201103, 0.3862607906, 0.3861016021, 0.3859425449, 0.3857836186, 0.3856248231, 0.3854661584, 0.3853076241, 0.3851492202, 0.3849909464, 0.3848328028, 0.3846747889, 0.3845169048, 0.3843591502, 0.3842015251, 0.3840440291, 0.3838866622, 0.3837294243, 0.3835723151, 0.3834153344, 0.3832584823, 0.3831017584, 0.3829451626, 0.3827886948, 0.3826323548, 0.3824761425, 0.3823200576, 0.3821641001, 0.3820082698, 0.3818525666, 0.3816969901, 0.3815415405, 0.3813862173, 0.3812310206, 0.3810759502, 0.3809258459, 0.3807710236, 0.380616327, 0.3804617562, 0.3803073108, 0.38015299071];
		var A2 = [55.19, 58.286, 58.238, 57.27, 55.556, 56.159, 55.317, 55.619, 55.492, 56.524, 55, 55.333, 54.19, 52.937, 51.317, 50.667, 52.143, 51.952, 51.683, 53, 52.73, 54.571, 53.222, 54.984, 55.444, 56.73, 58.54, 60.381, 59.73, 58.81, 55.683, 57.349, 58.016, 59.095, 59.857, 60.714, 59.571, 59.19, 56.444, 58.079, 56.921, 57.762, 57.73, 57.381, 56.444, 56.429, 58.27, 57.048, 55.937, 56.667, 58.524, 61.873, 61.683, 58.619, 57.651, 57.111, 57.968, 60.016, 57.365, 60.873, 62.524, 63.794, 62.714, 62.762, 63.333, 63.46, 62.206, 59.825, 59.46, 59.603, 58.27, 57.587, 57.222, 55.54, 58.952, 59.857, 59.714, 60.206, 59.794, 60.73, 59.746, 58.571, 58.222, 56.921, 56.476, 55.968, 55.984, 56.921, 58.175, 57.524, 57.556, 58.27, 58.952, 59.127, 61.46, 61.762, 59.937, 59.952, 59.794, 60.317, 61.429, 61.825, 59.508, 59.46, 59.27, 57.921, 57.587, 58.746, 57.571, 58.905, 60.651, 61.81, 62.381, 62.206, 60.444, 59.778, 60.111, 63.302, 64.365, 66.444, 67.238, 69.794, 69.603, 71.889, 72.429, 76.683, 78.46, 81.698, 83.27, 87.111, 85.746, 88.667, 90.571, 93.619, 95.984, 101.032, 104.222, 106.206, 111.238, 116.587, 118.429, 124.159, 128.762, 136.016, 141.286, 147.841, 153.714, 161.54, 168.286, 174.968, 182.651, 193.016, 200.079, 211, 221.556, 234.524, 246.81, 262.476, 280.222, 295.571, 311.032, 333.714, 353.905, 375.127, 398.333, 425.413, 453.222, 484.27, 516.333, 548.365, 583.968, 623.81, 666.825, 712.524, 756.952, 805, 857.048, 910.111, 965.968, 1024.619, 1087.444, 1153.984, 1221.635, 1295.857, 1375.762, 1451.952, 1535.413, 1625.397, 1716.254, 1809.127, 1910.54, 2013.984, 2124.333, 2238.73, 2362.175, 2491.524, 2619.81, 2751.778, 2893.921, 3041.524, 3197.508, 3349.937, 3509.841, 3674.349, 3838.317, 4008.302, 4192.381, 4378.111, 4561.397, 4749.873, 4942.143, 5142.619, 5346.159, 5548.476, 5760.381, 5968.492, 6176.032, 6397.492, 6623.651, 6847.873, 7068.683, 7285.063, 7513.73, 7739.587, 7961.778, 8199.778, 8443.937, 8680.286, 8907.619, 9144.762, 9386.143, 9624.111, 9859.048, 10096.032, 10328.698, 10558.857, 10788.873, 11024.222, 11260.063, 11491.286, 11718.778, 11950.841, 12176.111, 12399.619, 12621.857, 12848.317, 13056.27, 13239.968, 13434.889, 13643.825, 13840.413, 14028.206, 14231.286, 14434.857, 14636.159, 14842.032, 15023.603, 15209.317, 15393.127, 15570.397, 15723.841, 15884.365, 16027, 16184.032, 16347.587, 16502.032, 16644.508, 16781.333, 16931.984, 17047.937, 17152, 17269.889, 17375.619, 17479.206, 17565.175, 17635.698, 17692.683, 17759.905, 17827.889, 17861.714, 17893.698, 17955.889, 17978.429, 18032.746, 18068.349, 18049.286, 18052.81, 18046.063, 18072.254, 18067.73, 18077.127, 18089.921, 18100.159, 18116.46, 18130.127, 18132.651, 18126.27, 18156.444, 18182.587, 18199.032, 18250.508, 18297.175, 18324.143, 18376.46, 18398.302, 18435.73, 18499.127, 18554.968, 18622.095, 18672.048, 18741.063, 18820.873, 18875.603, 18911.302, 18972.937, 19014.111, 19066.746, 19140.889, 19198, 19259.381, 19311.937, 19366.587, 19413.46, 19432.222, 19456.254, 19479.587, 19501.19, 19547.794, 19571.714, 19581.19, 19570.714, 19589.095, 19558.984, 19571.079, 19572.381, 19551.81, 19534.825, 19515.222, 19496.619, 19458.063, 19415.413, 19380.778, 19344.635, 19305.841, 19299.952, 19274.81, 19253.921, 19215.127, 19146.143, 19102.508, 19065.651, 19043.508, 19001.619, 18984.063, 18927.73, 18898.825, 18861.095, 18825.571, 18792.889, 18763.079, 18717.032, 18683.921, 18634.063, 18556.46, 18500.794, 18423.333, 18380.19, 18326.571, 18266.429, 18236.048, 18179.587, 18117.937, 18051.857, 17981.032, 17914.349, 17860.238, 17792.413, 17726.143, 17663, 17600.238, 17529.175, 17487.048, 17420.127, 17369.302, 17340.349, 17294.714, 17277.127, 17236.841, 17201.413, 17196.524, 17140.714, 17109.619, 17118.921, 17107.27, 17074.19, 17076.825, 17074.54, 17086.968, 17100.968, 17105.683, 17117.397, 17140.73, 17138.683, 17164.73, 17183.921, 17219.286, 17269.873, 17310.19, 17365.54, 17429.016, 17477.111, 17542.46, 17602.825, 17628.032, 17673.889, 17740.603, 17799.683, 17869.984, 17921.492, 17969.46, 18016.841, 18066.46, 18104.016, 18148.667, 18189.698, 18226.016, 18224.206, 18204.365, 18196.825, 18155.556, 18119.905, 18091.571, 18061.587, 18039.222, 18011.365, 17950.206, 17876.444, 17790.937, 17705.016, 17621.19, 17550.032, 17484.968, 17421.175, 17359.111, 17288.683, 17204.063, 17107.937, 17022.841, 16921.667, 16826.873, 16759.635, 16679.46, 16599.857, 16518.794, 16432.825, 16347.238, 16294.968, 16238.603, 16173, 16105.365, 16030.238, 15939.048, 15834.365, 15732.254, 15628.524, 15529.032, 15429.81, 15360.841, 15303.619, 15280.206, 15224.333, 15171.365, 15102.984, 15008.698, 14914.254, 14836.508, 14762.349, 14672.095, 14606.27, 14555.984, 14494.81, 14436.063, 14393.317, 14353.778, 14312.175, 14249.587, 14219.032, 14174.143, 14121.016, 14052.524, 13963.571, 13895.492, 13824.286, 13786.397, 13754.048, 13725.079, 13658.159, 13594.429, 13556.571, 13499.921, 13414.81, 13343.873, 13291.381, 13211.556, 13153.46, 13055.524, 12983.349, 12887.524, 12805.381, 12708.111, 12605.619, 12520.302, 12428.476, 12312.619, 12188.524, 12069, 11963.349, 11852.365, 11720.683, 11576.175, 11443.016, 11308.762, 11177.302, 11017.19, 10874.794, 10722.444, 10595.365, 10443.905, 10293.159, 10145.698, 9992.905, 9847.952, 9687.079, 9514.937, 9361.286, 9196.381, 9039.175, 8871.746, 8724.159, 8598.873, 8456.095, 8308.46, 8180.841, 8029.048, 7881.762, 7727.397, 7576.921, 7433.032, 7298.619, 7157.54, 7012.175, 6876.794, 6754.587, 6630.27, 6509.444, 6383, 6251.683, 6113.921, 5988.143, 5850.54, 5712, 5586.222, 5461.111, 5328.762, 5207.365, 5094.873, 4995.841, 4885.746, 4769.333, 4657.032, 4549.032, 4450.063, 4357.079, 4248.476, 4142.286, 4047.984, 3953.825, 3868.683, 3785.095, 3694.048, 3610.476, 3532.175, 3454.968, 3375.238, 3294.032, 3218.048, 3143.397, 3065.937, 2990.698, 2919.73, 2851.476, 2787.302, 2717.762, 2646.81, 2577.54, 2512.603, 2442.571, 2376.889, 2311.667, 2250.286, 2194.603, 2136.651, 2079, 2019.984, 1964.381, 1914.587, 1861.349, 1812.937, 1760.857, 1708.016, 1656.873, 1606.921, 1555.667, 1507.222, 1462.508, 1418.746, 1375.857, 1333.143, 1290.524, 1248.762, 1204.619, 1166.111, 1127.492, 1084.683, 1046.587, 1010.413, 974.873, 941.571, 906.889, 876.063, 848.714, 819.222, 789.714, 761.159, 730.952, 702.937, 678.571, 654.778, 631.206, 611.19, 592.857, 572.889, 551, 532.365, 515.937, 499.667, 482.603, 469.476, 455.81, 442.429, 427.429, 414.127, 401.683, 389.619, 379.492, 369.825, 359.73, 352.889, 342.143, 333.619, 326.206, 317.429, 309.667, 303.762, 295.587, 290.54, 283.302, 276.873, 268.048, 263.302, 256.968, 252.286, 245.921, 240.778, 234.175, 230.587, 225.556, 220.413, 213.381, 208.587, 204.508, 199.206, 193.667, 192.603, 188.841, 186.048, 180.794, 176.73, 173.349, 171.365, 165.095, 163.794, 159.905, 155.937, 152.27, 150.238, 145.095, 143.952, 142.27, 140.254, 138.127, 135.27, 133.238, 130.762, 127.016, 122.206, 120.778, 119.127, 116.19, 114.825, 111.333, 111.063, 104.81, 102.27, 99.746, 97.111, 94.698, 94.444, 90.079, 87.254, 83.143, 82.794, 81.413, 76, 75.079, 69.857, 68.921, 64.889, 60.333, 59.873, 55.825, 53.492, 48.429, 47.683, 45.444, 42.857, 38.968, 37.619, 33.571, 33.095, 30.937, 29.968, 26.381, 23.048, 20.952, 18.587, 17.333, 17.032, 15.429, 13.159, 13.651, 12.81, 10.333, 9.873, 6.937, 7.556, 8.143, 10.254, 9.571, 9.619, 6.81, 6.079, 6.302, 6.571, 7.127, 9.206, 7.873, 11.46, 9.937, 11.905, 13.032, 11.333, 11.317, 12.921, 12.603, 14.937, 13.714, 14.19, 13.413, 14.444, 13.698, 15.016, 13.857, 16.048, 14.746, 14.762, 13.079, 13.556, 11.079, 12.794, 10, 13.571, 14.016, 14.016, 12.556, 15.397, 13.841, 14.46, 12.635, 12.111, 8.635, 7.825, 7.571, 7.524, 3.222, 3.444, 1.635, 2.984, 1.921, 3.238, 0.667, 2.841, -1.333, -3.905, -6.651, -6.286, -8.508, -7.27, -9.381, -7.667, -7.333, -3.238, -2.651, -3.079, -2.984, 1.19, 0.921, 3.286, 1.952, 3.111, 2.857, 3.556, 1.079, 2.238, 3.778, 6.143, 6.413, 7, 5.905, 4.063, 2.302, 2.508, -0.349, 0.095, 1.556, 2.651, -0.254, 0.778, -1.365, -0.286, 0.317, -1.698, -3.032, -1.873, -3.206, -5.27, -6.698, -7.429, -9.571, -10.381, -9.413, -9.952, -12.254, -12.111, -12.317, -13.683, -14.635, -14.984, -15.365, -13.683, -12.683, -12.905, -12.492, -9.492, -8.079, -6.524, -5.698, -5.016, -4.698, -3.127, -3.095, -3.159, -1.349, -2.143, -1.365, -2.175, -1.921, -1.921, -1.063, -0.302, -0.73, -1.032, -1.508, -2.349, -3.984, -4.857, -5.714, -6.714, -6.365, -4.286, -4.175, -1.413, -0.714, -1.238, 0.571, 0.873, 0.619, -0.746, -0.825, -0.635, -0.698, 0.492, -0.683, 0.492, 2.667, 4.143, 3.079, 2.952, 4.921, 5.921, 3.127, 4.302, 2.19, 3.349, 0.651, 0.54, -0.651, 1.508, 3.889, 4, 3.286, 2.19, 2.048, 3.143, 1.381, 0.857, -1.4137];
		bandlim.pulse = new bandlim.Pulse(
			f.reverse(),
			A2.reverse()
		);

	}

};
bandlim.init();
