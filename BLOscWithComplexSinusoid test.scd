// UGen BLOscWithComplexSinusoid only accepts control rate UGen as inputs.
// This UGen works, but CPU usage is much heavier than normal BLOsc. Therefore, it is not an efficient implementation (I made it for the research purpose).

// SynthDef
(
SynthDef(\blOscWithComplexSinusoidExample, { arg freq = 440.0, lo = 1, num = 15, slope = 1.0, eRatio = 1.0, amp = 1.0;
	Out.ar(0, BLOscWithComplexSinusoid.ar(freq, lo, num, slope, eRatio, mul: amp))
}).send
)

// Spectrum Analyzer
(
~openAnalyzer = {
var window, func;
window = Window("Analyzer", Rect(100, 100, 511, 300)); // width should be 511
func = FreqScopeView(window, window.view.bounds);
func.active_(true); // turn it on the first time;
func.freqMode_(1); //log scale
window.onClose_({ func.kill }); // you must have this
window.front;
}.value
)

// Demonstration with Gui
(
var window, sound, sl1, sl2, sl3, sl4, sl5, sl6;
window = Window("Band Limited Oscillator", Rect(100,420,620,300));
window.view.decorator = FlowLayout(window.view.bounds);
StaticText(window,620 @ 40).string_("Band Limited Oscillator").align_(\center)
	.font_(Font("Helvetica", 16));
window.view.decorator.nextLine;
sl1 = EZSlider(window, 600 @ 20, "Fundamental Freq", ControlSpec(55.0,880,\exponential,default: 440.0),
	 	{| view | if (sound.notNil)
	 			   { sound.set(\freq, view.value)}
	 	},
	labelWidth: 130);
window.view.decorator.nextLine;
sl2 = EZSlider(window, 600 @ 20, "Lowest Harmonics", ControlSpec(1,100,step: 1,default: 1),
	 	{| view | if (sound.notNil)
	 			   { sound.set(\lo, view.value)}
	 	},
	labelWidth: 130);
window.view.decorator.nextLine;
sl3 = EZSlider(window, 600 @ 20, "Num of Harmonics", ControlSpec(1,100,step: 1,default: 7),
	 	{| view | if (sound.notNil)
	 			   { sound.set(\num, view.value)}
	 	},
	labelWidth: 130);
window.view.decorator.nextLine;
sl4 = EZSlider(window, 600 @ 20, "Slope", ControlSpec(0.1,2.0,default: 1.0),
	 	{| view | if (sound.notNil)
	 			   { sound.set(\slope, view.value)}
	 	},
	labelWidth: 130);
window.view.decorator.nextLine;
sl5 = EZSlider(window, 600 @ 20, "Even to Odd Ratio", ControlSpec(0,1.0,default: 1.0),
	 	{| view | if (sound.notNil)
	 			   { sound.set(\eRatio, view.value)}
	 	},
	labelWidth: 130);
window.view.decorator.nextLine;
sl6 = EZSlider(window, 600 @ 20, "Amp", ControlSpec(0,1.0,\amp,default: 1.0),
	 	{| view | if (sound.notNil)
	 			   { sound.set(\amp, view.value)}
	 	},
	labelWidth: 130);
window.view.decorator.left_(300).top_(230);
Button(window, 50@50)
	.states_([["Start"], ["Stop", Color.black, Color.red(alpha: 0.2)]])
	.action_({|view|  if (view.value ==1 && sound.isNil,
	{sound =
		Synth(\blOscWithComplexSinusoidExample,
			[\freq,sl1.value,
			\lo, sl2.value,
			\num, sl3.value,
			\slope, sl4.value,
			\eRatio, sl5.value,
			\amp, sl6.value])},
						{sound.free; sound = nil }
						)
		});
window.front;
)