/*
 Copyright (c) 2021 Temporal Guild Group, Austral University of Chile, Valdivia Chile.
 This file and all powermon software is licensed under the MIT License. 
 Please refer to LICENSE for more details.
 */
#include <unistd.h>
#include <cstdint>
#include <cstring>

#ifndef RAPL_H_
#define RAPL_H_

#define MAX_LINE 256
#define MAX_CPU 256
#define MAX_SOCKETS 8

#define COOLDOWN_MS  20
struct rapl_state_t {
	uint64_t pkg;
	uint64_t pp0;
	uint64_t pp1;
	uint64_t dram;
	struct timeval tsc;
};

// CPU power measure functions
void CPUPowerBegin(const char *alg, int ms);
void CPUPowerEnd();

// pthread functions
void *CPUpowerPollingFunc(void *ptr);

class Rapl {

private:
	// Rapl configuration
	int fd[MAX_SOCKETS];
	int core = 0;
	bool pp1_supported = true;
	//vendor 0=Intel, 1=AMD
	int vendor;
	int n_sockets;
	int smt;
	int n_logical_cores;
	int sockets[MAX_SOCKETS];
	int first_lcoreid[MAX_SOCKETS];
	double power_units, energy_units, time_units;
	double thermal_spec_power, minimum_power, maximum_power, time_window;

	// Rapl state
	rapl_state_t *current_state[MAX_SOCKETS];
	rapl_state_t *prev_state[MAX_SOCKETS];
	rapl_state_t *next_state[MAX_SOCKETS];
	rapl_state_t *running_total[MAX_SOCKETS];

	bool detect_pp1();
	int get_vendor();
	void open_msr(int socket, int core);
	uint64_t read_msr(int socket, uint32_t msr_offset);
	double time_delta(struct timeval *begin, struct timeval *after);
	uint64_t energy_delta(uint64_t before, uint64_t after);
	double power(uint64_t before, uint64_t after, double time_delta);

public:
	Rapl();
	void reset();
	void sample();
	void sample(int socket);

	double pkg_current_power();
	double pp0_current_power();
	double pp1_current_power();
	double dram_current_power();

	double pkg_average_power();
	double pp0_average_power();
	double pp1_average_power();
	double dram_average_power();

	double pkg_total_energy();
	double pp0_total_energy();
	double pp1_total_energy();
	double dram_total_energy();

	double total_time();
	double current_time();
	int get_n_sockets();
	int get_n_logical_cores();
	int get_smt();
};

#endif /* RAPL_H_ */
