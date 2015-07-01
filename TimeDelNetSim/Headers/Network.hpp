#ifndef NETWORK_H
#define NETWORK_H
struct Synapse{
	int		NStart;
	int		NEnd;
	float	Weight;
	int	DelayinTsteps;
};

struct Neuron{
	float a;
	float b;
	float c;
	float d;
};
#endif