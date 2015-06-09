//
//  easing.c
//
//  Copyright (c) 2011, Auerhaus Development, LLC
//
//  This program is free software. It comes without any warranty, to
//  the extent permitted by applicable law. You can redistribute it
//  and/or modify it under the terms of the Do What The Fuck You Want
//  To Public License, Version 2, as published by Sam Hocevar. See
//  http://sam.zoy.org/wtfpl/COPYING for more details.
//

#include <math.h>

float backEaseInOutf(float p);
float backEaseInf(float p);
float backEaseOutf(float p);
float bounceEaseInOutf(float p);
float bounceEaseInf(float p);
float bounceEaseOutf(float p);
float circularEaseInOutf(float p);
float circularEaseInf(float p);
float circularEaseOutf(float p);
float cubicEaseInOutf(float p);
float cubicEaseInf(float p);
float cubicEaseOutf(float p);
float elasticEaseInOutf(float p);
float elasticEaseInf(float p);
float elasticEaseOutf(float p);
float exponentialEaseInOutf(float p);
float exponentialEaseInf(float p);
float exponentialEaseOutf(float p);
float linearInterpolationf(float p);
float quadraticEaseInOutf(float p);
float quadraticEaseInf(float p);
float quadraticEaseOutf(float p);
float quarticEaseInOutf(float p);
float quarticEaseInf(float p);
float quarticEaseOutf(float p);
float quinticEaseInOutf(float p);
float quinticEaseInf(float p);
float quinticEaseOutf(float p);
float sineEaseInOutf(float p);
float sineEaseInf(float p);
float sineEaseOutf(float p);
double backEaseInOutd(double p);
double backEaseInd(double p);
double backEaseOutd(double p);
double bounceEaseInOutd(double p);
double bounceEaseInd(double p);
double bounceEaseOutd(double p);
double circularEaseInOutd(double p);
double circularEaseInd(double p);
double circularEaseOutd(double p);
double cubicEaseInOutd(double p);
double cubicEaseInd(double p);
double cubicEaseOutd(double p);
double elasticEaseInOutd(double p);
double elasticEaseInd(double p);
double elasticEaseOutd(double p);
double exponentialEaseInOutd(double p);
double exponentialEaseInd(double p);
double exponentialEaseOutd(double p);
double linearInterpolationd(double p);
double quadraticEaseInOutd(double p);
double quadraticEaseInd(double p);
double quadraticEaseOutd(double p);
double quarticEaseInOutd(double p);
double quarticEaseInd(double p);
double quarticEaseOutd(double p);
double quinticEaseInOutd(double p);
double quinticEaseInd(double p);
double quinticEaseOutd(double p);
double sineEaseInOutd(double p);
double sineEaseInd(double p);
double sineEaseOutd(double p);

// Modeled after the line y = x
float linearInterpolationf(float p)
{
	return p;
}

// Modeled after the parabola y = x^2
float quadraticEaseInf(float p)
{
	return p * p;
}

// Modeled after the parabola y = -x^2 + 2x
float quadraticEaseOutf(float p)
{
	return -(p * (p - 2));
}

// Modeled after the piecewise quadratic
// y = (1/2)((2x)^2)             ; [0, 0.5)
// y = -(1/2)((2x-1)*(2x-3) - 1) ; [0.5, 1]
float quadraticEaseInOutf(float p)
{
	if(p < 0.5)
	{
		return 2 * p * p;
	}
	else
	{
		return (-2 * p * p) + (4 * p) - 1;
	}
}

// Modeled after the cubic y = x^3
float cubicEaseInf(float p)
{
	return p * p * p;
}

// Modeled after the cubic y = (x - 1)^3 + 1
float cubicEaseOutf(float p)
{
	float f = (p - 1);
	return f * f * f + 1;
}

// Modeled after the piecewise cubic
// y = (1/2)((2x)^3)       ; [0, 0.5)
// y = (1/2)((2x-2)^3 + 2) ; [0.5, 1]
float cubicEaseInOutf(float p)
{
	if(p < 0.5)
	{
		return 4 * p * p * p;
	}
	else
	{
		float f = ((2 * p) - 2);
		return 0.5 * f * f * f + 1;
	}
}

// Modeled after the quartic x^4
float quarticEaseInf(float p)
{
	return p * p * p * p;
}

// Modeled after the quartic y = 1 - (x - 1)^4
float quarticEaseOutf(float p)
{
	float f = (p - 1);
	return f * f * f * (1 - p) + 1;
}

// Modeled after the piecewise quartic
// y = (1/2)((2x)^4)        ; [0, 0.5)
// y = -(1/2)((2x-2)^4 - 2) ; [0.5, 1]
float quarticEaseInOutf(float p) 
{
	if(p < 0.5)
	{
		return 8 * p * p * p * p;
	}
	else
	{
		float f = (p - 1);
		return -8 * f * f * f * f + 1;
	}
}

// Modeled after the quintic y = x^5
float quinticEaseInf(float p) 
{
	return p * p * p * p * p;
}

// Modeled after the quintic y = (x - 1)^5 + 1
float quinticEaseOutf(float p) 
{
	float f = (p - 1);
	return f * f * f * f * f + 1;
}

// Modeled after the piecewise quintic
// y = (1/2)((2x)^5)       ; [0, 0.5)
// y = (1/2)((2x-2)^5 + 2) ; [0.5, 1]
float quinticEaseInOutf(float p) 
{
	if(p < 0.5)
	{
		return 16 * p * p * p * p * p;
	}
	else
	{
		float f = ((2 * p) - 2);
		return  0.5 * f * f * f * f * f + 1;
	}
}

// Modeled after quarter-cycle of sine wave
float sineEaseInf(float p)
{
	return sin((p - 1) * M_PI_2) + 1;
}

// Modeled after quarter-cycle of sine wave (different phase)
float sineEaseOutf(float p)
{
	return sin(p * M_PI_2);
}

// Modeled after half sine wave
float sineEaseInOutf(float p)
{
	return 0.5 * (1 - cos(p * M_PI));
}

// Modeled after shifted quadrant IV of unit circle
float circularEaseInf(float p)
{
	return 1 - sqrt(1 - (p * p));
}

// Modeled after shifted quadrant II of unit circle
float circularEaseOutf(float p)
{
	return sqrt((2 - p) * p);
}

// Modeled after the piecewise circular function
// y = (1/2)(1 - sqrt(1 - 4x^2))           ; [0, 0.5)
// y = (1/2)(sqrt(-(2x - 3)*(2x - 1)) + 1) ; [0.5, 1]
float circularEaseInOutf(float p)
{
	if(p < 0.5)
	{
		return 0.5 * (1 - sqrt(1 - 4 * (p * p)));
	}
	else
	{
		return 0.5 * (sqrt(-((2 * p) - 3) * ((2 * p) - 1)) + 1);
	}
}

// Modeled after the exponential function y = 2^(10(x - 1))
float exponentialEaseInf(float p)
{
	return (p == 0.0) ? p : pow(2, 10 * (p - 1));
}

// Modeled after the exponential function y = -2^(-10x) + 1
float exponentialEaseOutf(float p)
{
	return (p == 1.0) ? p : 1 - pow(2, -10 * p);
}

// Modeled after the piecewise exponential
// y = (1/2)2^(10(2x - 1))         ; [0,0.5)
// y = -(1/2)*2^(-10(2x - 1))) + 1 ; [0.5,1]
float exponentialEaseInOutf(float p)
{
	if(p == 0.0 || p == 1.0) return p;
	
	if(p < 0.5)
	{
		return 0.5 * pow(2, (20 * p) - 10);
	}
	else
	{
		return -0.5 * pow(2, (-20 * p) + 10) + 1;
	}
}

// Modeled after the damped sine wave y = sin(13pi/2*x)*pow(2, 10 * (x - 1))
float elasticEaseInf(float p)
{
	return sin(13 * M_PI_2 * p) * pow(2, 10 * (p - 1));
}

// Modeled after the damped sine wave y = sin(-13pi/2*(x + 1))*pow(2, -10x) + 1
float elasticEaseOutf(float p)
{
	return sin(-13 * M_PI_2 * (p + 1)) * pow(2, -10 * p) + 1;
}

// Modeled after the piecewise exponentially-damped sine wave:
// y = (1/2)*sin(13pi/2*(2*x))*pow(2, 10 * ((2*x) - 1))      ; [0,0.5)
// y = (1/2)*(sin(-13pi/2*((2x-1)+1))*pow(2,-10(2*x-1)) + 2) ; [0.5, 1]
float elasticEaseInOutf(float p)
{
	if(p < 0.5)
	{
		return 0.5 * sin(13 * M_PI_2 * (2 * p)) * pow(2, 10 * ((2 * p) - 1));
	}
	else
	{
		return 0.5 * (sin(-13 * M_PI_2 * ((2 * p - 1) + 1)) * pow(2, -10 * (2 * p - 1)) + 2);
	}
}

// Modeled after the overshooting cubic y = x^3-x*sin(x*pi)
float backEaseInf(float p)
{
	return p * p * p - p * sin(p * M_PI);
}

// Modeled after overshooting cubic y = 1-((1-x)^3-(1-x)*sin((1-x)*pi))
float backEaseOutf(float p)
{
	float f = (1 - p);
	return 1 - (f * f * f - f * sin(f * M_PI));
}

// Modeled after the piecewise overshooting cubic function:
// y = (1/2)*((2x)^3-(2x)*sin(2*x*pi))           ; [0, 0.5)
// y = (1/2)*(1-((1-x)^3-(1-x)*sin((1-x)*pi))+1) ; [0.5, 1]
float backEaseInOutf(float p)
{
	if(p < 0.5)
	{
		float f = 2 * p;
		return 0.5 * (f * f * f - f * sin(f * M_PI));
	}
	else
	{
		float f = (1 - (2*p - 1));
		return 0.5 * (1 - (f * f * f - f * sin(f * M_PI))) + 0.5;
	}
}

float bounceEaseInf(float p)
{
	return 1 - bounceEaseOutf(1 - p);
}

float bounceEaseOutf(float p)
{
	if(p < 4/11.0)
	{
		return (121 * p * p)/16.0;
	}
	else if(p < 8/11.0)
	{
		return (363/40.0 * p * p) - (99/10.0 * p) + 17/5.0;
	}
	else if(p < 9/10.0)
	{
		return (4356/361.0 * p * p) - (35442/1805.0 * p) + 16061/1805.0;
	}
	else
	{
		return (54/5.0 * p * p) - (513/25.0 * p) + 268/25.0;
	}
}

float bounceEaseInOutf(float p)
{
	if(p < 0.5)
	{
		return 0.5 * bounceEaseInf(p*2);
	}
	else
	{
		return 0.5 * bounceEaseOutf(p * 2 - 1) + 0.5;
	}
}



// Modeled after the line y = x
double linearInterpolationd(double p)
{
	return p;
}

// Modeled after the parabola y = x^2
double quadraticEaseInd(double p)
{
	return p * p;
}

// Modeled after the parabola y = -x^2 + 2x
double quadraticEaseOutd(double p)
{
	return -(p * (p - 2));
}

// Modeled after the piecewise quadratic
// y = (1/2)((2x)^2)             ; [0, 0.5)
// y = -(1/2)((2x-1)*(2x-3) - 1) ; [0.5, 1]
double quadraticEaseInOutd(double p)
{
	if(p < 0.5)
	{
		return 2 * p * p;
	}
	else
	{
		return (-2 * p * p) + (4 * p) - 1;
	}
}

// Modeled after the cubic y = x^3
double cubicEaseInd(double p)
{
	return p * p * p;
}

// Modeled after the cubic y = (x - 1)^3 + 1
double cubicEaseOutd(double p)
{
	double f = (p - 1);
	return f * f * f + 1;
}

// Modeled after the piecewise cubic
// y = (1/2)((2x)^3)       ; [0, 0.5)
// y = (1/2)((2x-2)^3 + 2) ; [0.5, 1]
double cubicEaseInOutd(double p)
{
	if(p < 0.5)
	{
		return 4 * p * p * p;
	}
	else
	{
		double f = ((2 * p) - 2);
		return 0.5 * f * f * f + 1;
	}
}

// Modeled after the quartic x^4
double quarticEaseInd(double p)
{
	return p * p * p * p;
}

// Modeled after the quartic y = 1 - (x - 1)^4
double quarticEaseOutd(double p)
{
	double f = (p - 1);
	return f * f * f * (1 - p) + 1;
}

// Modeled after the piecewise quartic
// y = (1/2)((2x)^4)        ; [0, 0.5)
// y = -(1/2)((2x-2)^4 - 2) ; [0.5, 1]
double quarticEaseInOutd(double p) 
{
	if(p < 0.5)
	{
		return 8 * p * p * p * p;
	}
	else
	{
		double f = (p - 1);
		return -8 * f * f * f * f + 1;
	}
}

// Modeled after the quintic y = x^5
double quinticEaseInd(double p) 
{
	return p * p * p * p * p;
}

// Modeled after the quintic y = (x - 1)^5 + 1
double quinticEaseOutd(double p) 
{
	double f = (p - 1);
	return f * f * f * f * f + 1;
}

// Modeled after the piecewise quintic
// y = (1/2)((2x)^5)       ; [0, 0.5)
// y = (1/2)((2x-2)^5 + 2) ; [0.5, 1]
double quinticEaseInOutd(double p) 
{
	if(p < 0.5)
	{
		return 16 * p * p * p * p * p;
	}
	else
	{
		double f = ((2 * p) - 2);
		return  0.5 * f * f * f * f * f + 1;
	}
}

// Modeled after quarter-cycle of sine wave
double sineEaseInd(double p)
{
	return sin((p - 1) * M_PI_2) + 1;
}

// Modeled after quarter-cycle of sine wave (different phase)
double sineEaseOutd(double p)
{
	return sin(p * M_PI_2);
}

// Modeled after half sine wave
double sineEaseInOutd(double p)
{
	return 0.5 * (1 - cos(p * M_PI));
}

// Modeled after shifted quadrant IV of unit circle
double circularEaseInd(double p)
{
	return 1 - sqrt(1 - (p * p));
}

// Modeled after shifted quadrant II of unit circle
double circularEaseOutd(double p)
{
	return sqrt((2 - p) * p);
}

// Modeled after the piecewise circular function
// y = (1/2)(1 - sqrt(1 - 4x^2))           ; [0, 0.5)
// y = (1/2)(sqrt(-(2x - 3)*(2x - 1)) + 1) ; [0.5, 1]
double circularEaseInOutd(double p)
{
	if(p < 0.5)
	{
		return 0.5 * (1 - sqrt(1 - 4 * (p * p)));
	}
	else
	{
		return 0.5 * (sqrt(-((2 * p) - 3) * ((2 * p) - 1)) + 1);
	}
}

// Modeled after the exponential function y = 2^(10(x - 1))
double exponentialEaseInd(double p)
{
	return (p == 0.0) ? p : pow(2, 10 * (p - 1));
}

// Modeled after the exponential function y = -2^(-10x) + 1
double exponentialEaseOutd(double p)
{
	return (p == 1.0) ? p : 1 - pow(2, -10 * p);
}

// Modeled after the piecewise exponential
// y = (1/2)2^(10(2x - 1))         ; [0,0.5)
// y = -(1/2)*2^(-10(2x - 1))) + 1 ; [0.5,1]
double exponentialEaseInOutd(double p)
{
	if(p == 0.0 || p == 1.0) return p;
	
	if(p < 0.5)
	{
		return 0.5 * pow(2, (20 * p) - 10);
	}
	else
	{
		return -0.5 * pow(2, (-20 * p) + 10) + 1;
	}
}

// Modeled after the damped sine wave y = sin(13pi/2*x)*pow(2, 10 * (x - 1))
double elasticEaseInd(double p)
{
	return sin(13 * M_PI_2 * p) * pow(2, 10 * (p - 1));
}

// Modeled after the damped sine wave y = sin(-13pi/2*(x + 1))*pow(2, -10x) + 1
double elasticEaseOutd(double p)
{
	return sin(-13 * M_PI_2 * (p + 1)) * pow(2, -10 * p) + 1;
}

// Modeled after the piecewise exponentially-damped sine wave:
// y = (1/2)*sin(13pi/2*(2*x))*pow(2, 10 * ((2*x) - 1))      ; [0,0.5)
// y = (1/2)*(sin(-13pi/2*((2x-1)+1))*pow(2,-10(2*x-1)) + 2) ; [0.5, 1]
double elasticEaseInOutd(double p)
{
	if(p < 0.5)
	{
		return 0.5 * sin(13 * M_PI_2 * (2 * p)) * pow(2, 10 * ((2 * p) - 1));
	}
	else
	{
		return 0.5 * (sin(-13 * M_PI_2 * ((2 * p - 1) + 1)) * pow(2, -10 * (2 * p - 1)) + 2);
	}
}

// Modeled after the overshooting cubic y = x^3-x*sin(x*pi)
double backEaseInd(double p)
{
	return p * p * p - p * sin(p * M_PI);
}

// Modeled after overshooting cubic y = 1-((1-x)^3-(1-x)*sin((1-x)*pi))
double backEaseOutd(double p)
{
	double f = (1 - p);
	return 1 - (f * f * f - f * sin(f * M_PI));
}

// Modeled after the piecewise overshooting cubic function:
// y = (1/2)*((2x)^3-(2x)*sin(2*x*pi))           ; [0, 0.5)
// y = (1/2)*(1-((1-x)^3-(1-x)*sin((1-x)*pi))+1) ; [0.5, 1]
double backEaseInOutd(double p)
{
	if(p < 0.5)
	{
		double f = 2 * p;
		return 0.5 * (f * f * f - f * sin(f * M_PI));
	}
	else
	{
		double f = (1 - (2*p - 1));
		return 0.5 * (1 - (f * f * f - f * sin(f * M_PI))) + 0.5;
	}
}

double bounceEaseInd(double p)
{
	return 1 - bounceEaseOutd(1 - p);
}

double bounceEaseOutd(double p)
{
	if(p < 4/11.0)
	{
		return (121 * p * p)/16.0;
	}
	else if(p < 8/11.0)
	{
		return (363/40.0 * p * p) - (99/10.0 * p) + 17/5.0;
	}
	else if(p < 9/10.0)
	{
		return (4356/361.0 * p * p) - (35442/1805.0 * p) + 16061/1805.0;
	}
	else
	{
		return (54/5.0 * p * p) - (513/25.0 * p) + 268/25.0;
	}
}

double bounceEaseInOutd(double p)
{
	if(p < 0.5)
	{
		return 0.5 * bounceEaseInd(p*2);
	}
	else
	{
		return 0.5 * bounceEaseOutd(p * 2 - 1) + 0.5;
	}
}

