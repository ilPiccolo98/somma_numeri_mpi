#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#define SHIFT_FOR_ARGV 2
#define INPUT_NUMBERS_LIMIT 20
#define OFFSET_STRATEGY 2
#define OFFSET_ID_PRINTER_PROCESS 1
#define FIXED_PARAMETERS 4

int get_size_numbers_to_sum(char **argv)
{
	int size = atoi(argv[1]);
	return size;
}

float* get_numbers_from_argv(int size, char *argv[])
{
	float *numbers_to_sum = (float*)calloc(size, sizeof(float));
	int i;
	for(i = 0; i != size; ++i)
		numbers_to_sum[i] = atof(argv[i + SHIFT_FOR_ARGV]);
	return numbers_to_sum;
}

float* get_random_numbers(int size)
{
	float *numbers_to_sum = (float*)calloc(size, sizeof(float));
	int i;
	for(i = 0; i != size; ++i)
		numbers_to_sum[i] = rand();
	return numbers_to_sum;
}

float* get_numbers_to_sum(int size, char *argv[])
{
	return (size <= INPUT_NUMBERS_LIMIT) ? get_numbers_from_argv(size, argv) : get_random_numbers(size);
}

int get_strategy(char *argv[], int argc)
{
	return atoi(argv[argc - OFFSET_STRATEGY]);
}

int get_id_printer_processor(char *argv[], int argc)
{
	return atoi(argv[argc - OFFSET_ID_PRINTER_PROCESS]);
}

void send_numbers(float *local_numbers_to_sum, int local_size_numbers_to_sum, int processors, int rest)
{
	int id_comunication = 22;
	int step = local_size_numbers_to_sum;
	int start = 0;
	int i;
	for(i = 1; i < processors; ++i)
	{
		start += step;
		id_comunication = 22 + i;
		if(i == rest)
			step -= 1;
		MPI_Send(&local_numbers_to_sum[start], step, MPI_FLOAT, i, id_comunication, MPI_COMM_WORLD);
	}
}

void get_numbers(float *local_numbers_to_sum, int local_size_numbers_to_sum, int id_processor)
{
	int id_comunication = 22 + id_processor;
	MPI_Status status;
	MPI_Recv(local_numbers_to_sum, local_size_numbers_to_sum, MPI_FLOAT, 0, id_comunication, MPI_COMM_WORLD, &status);
}

float get_sum(float *local_numbers_to_sum, int local_size_numbers_to_sum)
{
	int i;
	float sum = 0;
	for(i = 0; i < local_size_numbers_to_sum; ++i)
		sum += local_numbers_to_sum[i];
	return sum;
}

int is_a_power_of_2(int number)
{
	return number > 1 && ((number & (number-1)) == 0);
}

int get_tree_hight(int processors)
{
	return (int)(log10(processors)/log10(2));
}

void strategy1(float *sum, int id_processor, int processors, double t0)
{
	if (id_processor == 0)
	{
		int i;
		for(i = 1; i < processors; ++i)
		{
			MPI_Status status;
			int id_comunication = 80 + i;
			float partial_sum = 0;
			MPI_Recv(&partial_sum, 1, MPI_FLOAT, i, id_comunication, MPI_COMM_WORLD, &status);
			double end = MPI_Wtime();
			*sum += partial_sum;
		}
	}
	else
	{
		int id_comunication = id_processor + 80;
		MPI_Send(sum, 1, MPI_FLOAT, 0, id_comunication, MPI_COMM_WORLD);
	}
}

void strategy2(float *sum, int id_processor, int processors)
{
	int i;
	for(i = 0; i < get_tree_hight(processors); ++i)
	{
		if(id_processor % (int)pow(2, i) == 0)
		{
			if(id_processor % (int)pow(2, i + 1) == 0)
			{
				int id_comunication = 100 + i*10 + id_processor;
				float partial_sum = 0;
				MPI_Status status;
				MPI_Recv(&partial_sum, 1, MPI_FLOAT, id_processor + (int)pow(2, i), id_comunication, MPI_COMM_WORLD, &status);
				*sum += partial_sum;				
			}
			else
			{	
				int id_comunication = 100 + i*10 + id_processor - (int)pow(2, i);
				MPI_Send(sum, 1, MPI_FLOAT, id_processor - (int)pow(2, i), id_comunication, MPI_COMM_WORLD);
			}
		}
	}
}

void strategy3(float *sum, int id_processor, int processors)
{
	int i;
	for(i = 0; i < get_tree_hight(processors); ++i)
	{
		if(id_processor % (int)pow(2, i + 1) < (int)pow(2, i))
		{
			MPI_Status status;
			float received_partial_sum = 0;
			int id_comunication = 200 + i * 10 + id_processor;
			MPI_Recv(&received_partial_sum, 1, MPI_FLOAT, id_processor + (int)pow(2, i), id_comunication, MPI_COMM_WORLD, &status);
			id_comunication = 200 + id_processor + (int)pow(2, i) + i * 10;
			MPI_Send(sum, 1, MPI_FLOAT, id_processor + (int)pow(2, i), id_comunication, MPI_COMM_WORLD);
			*sum += received_partial_sum;
		}
		else
		{
			MPI_Status status;
			float received_partial_sum = 0;
			int id_comunication = 200 + i * 10 + id_processor - (int)pow(2, i);
			MPI_Send(sum, 1, MPI_FLOAT, id_processor - (int)pow(2, i), id_comunication, MPI_COMM_WORLD);
			id_comunication = 200 + i * 10 + id_processor;
			MPI_Recv(&received_partial_sum, 1, MPI_FLOAT, id_processor - (int)pow(2, i), id_comunication, MPI_COMM_WORLD, &status);
			*sum += received_partial_sum;
		}
	}
}

void sum_numbers(int strategy, float *sum, int size_numbers_to_sum, int id_processor, int processors, double t0)
{
	if(strategy == 1 || !is_a_power_of_2(processors))
	{
		if(id_processor == 0 && strategy != 1)
			puts("I'm using the strategy 1. The number of processors is not a power of 2");
		strategy1(sum, id_processor, processors, t0);
	}
	else if(strategy == 2)
		strategy2(sum, id_processor, processors);
	else if(strategy == 3)
		strategy3(sum, id_processor, processors);
}

void free_memory(int id_processor, float *numbers_to_sum, float *local_numbers_to_sum)
{
	if(id_processor == 0)
		free(numbers_to_sum);
	else
		free(local_numbers_to_sum);
}

void print_sum(int id_processor, int id_printer_processor, float sum, double total_time_elapsed)
{
	if(id_printer_processor == -1 || id_printer_processor == id_processor)
		printf("I'm the processor %d - This is the sum %f - Total time elapsed: %f seconds\n", id_processor, sum, total_time_elapsed);
}

void send_sum(int id_printer_processor, int id_processor, float *sum)
{
	if(id_printer_processor != -1 && id_printer_processor != 0)
	{
		if(id_processor == 0)
		{
			int id_comunication = 300 + id_printer_processor;
			MPI_Send(sum, 1, MPI_FLOAT, id_printer_processor, id_comunication, MPI_COMM_WORLD);
		}
		else if(id_processor == id_printer_processor)
		{
			MPI_Status status;
			int id_comunication = 300 + id_printer_processor;
			MPI_Recv(sum, 1, MPI_FLOAT, 0, id_comunication, MPI_COMM_WORLD, &status);
		}
	}
	else if(id_printer_processor == -1)
		MPI_Bcast(sum, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
}

void send_total_time_elapsed(int id_printer_processor, int id_processor, double *total_time_elapsed)
{
	if(id_printer_processor != -1 && id_printer_processor != 0)
	{
		if(id_processor == 0)
		{
			int id_comunication = 400 + id_printer_processor;
			MPI_Send(total_time_elapsed, 1, MPI_DOUBLE, id_printer_processor, id_comunication, MPI_COMM_WORLD);
		}
		else if(id_processor == id_printer_processor)
		{
			MPI_Status status;
			int id_comunication = 400 + id_printer_processor;
			MPI_Recv(total_time_elapsed, 1, MPI_DOUBLE, 0, id_comunication, MPI_COMM_WORLD, &status);
		}
	}
	else if(id_printer_processor == -1)
		MPI_Bcast(total_time_elapsed, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void handle_error(char *error_message)
{
	fprintf(stderr, error_message);
	MPI_Abort(MPI_COMM_WORLD, -1);
	exit(-1);
}

int main(int argc, char *argv[])
{
	srand((unsigned int)time(NULL));
	MPI_Init(&argc, &argv);
	int id_processor = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &id_processor);
	int processors = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &processors);
	int size_numbers_to_sum = 0;
	int local_size_numbers_to_sum = 0;
	float *numbers_to_sum = NULL;
	float *local_numbers_to_sum = NULL;
	int strategy = 1;
	int id_printer_processor = -1;
	int id_comunication = 0;
	if(id_processor == 0)
	{
		if(argc < FIXED_PARAMETERS)
		{
			free_memory(id_processor, numbers_to_sum, local_numbers_to_sum);
			handle_error("Error number of arguments!");
		}
		size_numbers_to_sum = get_size_numbers_to_sum(argv);
		if(size_numbers_to_sum > 1 || (size_numbers_to_sum <= INPUT_NUMBERS_LIMIT && argc == size_numbers_to_sum + FIXED_PARAMETERS) || (size_numbers_to_sum > INPUT_NUMBERS_LIMIT && argc == FIXED_PARAMETERS))
			numbers_to_sum = get_numbers_to_sum(size_numbers_to_sum, argv);
		else
		{
			free_memory(id_processor, numbers_to_sum, local_numbers_to_sum);
			handle_error("Error size numbers input!");
		}
		strategy = get_strategy(argv, argc);
		if(strategy < 1 || strategy > 3)
		{
			free_memory(id_processor, numbers_to_sum, local_numbers_to_sum);
			handle_error("Error strategy!");
		}
		id_printer_processor = get_id_printer_processor(argv, argc);
		if(id_printer_processor < -1 || id_printer_processor >= processors)
		{
			free_memory(id_processor, numbers_to_sum, local_numbers_to_sum);
			handle_error("Error id printer processor!");		
		}
	}
	MPI_Bcast(&size_numbers_to_sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	local_size_numbers_to_sum = size_numbers_to_sum / processors;
	int rest = size_numbers_to_sum % processors;
	if (id_processor < rest) 
	    ++local_size_numbers_to_sum;
	if(id_processor == 0)
	{
		local_numbers_to_sum = numbers_to_sum;
		send_numbers(local_numbers_to_sum, local_size_numbers_to_sum, processors, rest);
	}
	else 
	{
		local_numbers_to_sum = (float*)calloc(local_size_numbers_to_sum, sizeof(float));
		get_numbers(local_numbers_to_sum, local_size_numbers_to_sum, id_processor);
	}
	MPI_Bcast(&strategy, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	double t0 = MPI_Wtime();
	float sum = get_sum(local_numbers_to_sum, local_size_numbers_to_sum);
	sum_numbers(strategy, &sum, size_numbers_to_sum, id_processor, processors, t0);
	double t1 = MPI_Wtime();
	double time = t1 - t0;
	MPI_Bcast(&id_printer_processor, 1, MPI_INT, 0, MPI_COMM_WORLD);
	send_sum(id_printer_processor, id_processor, &sum);
	double total_time_elapsed = 0;
	MPI_Reduce(&time, &total_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	send_total_time_elapsed(id_printer_processor, id_processor, &total_time_elapsed);
	print_sum(id_processor, id_printer_processor, sum, total_time_elapsed);
	free_memory(id_processor, numbers_to_sum, local_numbers_to_sum);
	MPI_Finalize();
	return 0;
}
