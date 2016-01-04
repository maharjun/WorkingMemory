#ifndef NETWORK_H
#define NETWORK_H
#include <cstdint>
struct Synapse{
	int32_t NStart;
	int32_t NEnd;
	float	Weight;
	int32_t	DelayinTsteps;
};

struct Neuron{
	float a;
	float b;
	float c;
	float d;
};
#endif