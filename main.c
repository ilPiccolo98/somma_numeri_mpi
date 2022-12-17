// ----------------------------------------------------------------------------------------------------
// | SCOPO                                                                                            
// | PROGRAMMA PER ESEGUIRE LA SOMMA DEI NUMERI IN PARALLELO SU UN'ARCHITETTURA A MEMORIA DISTRIBUITA 
// | TRAMITE LA LIBRERIA MPI. IL PROGRAMMA PUO' SVOLGERE LA SOMMA TRAMITE 3 DIVERSE STRATEGIE. LA     
// | STRATEGIA DEVE ESSERE SPECIFICATA QUANDO VIENE ESEGUITO IL COMANDO mpiexec COME PARAMETRO.       
// | I NUMERI DA SOMMARE DEVONO ESSERE SPECIFICATI SE E SOLO SE LA QUANTITA' DEI NUMERI DICHIARATA    
// | E' <= 20, SEMPRE TRAMITE PARAMETRO. PUO' ESSERE SPECIFICATO IL PROCESSORE CHE DEVE STAMPARE      
// | LA SOMMA FINALE E IL TEMPO DI ESECUZIONE, OPPURE SI PUO' SPECIFICARE CHE TUTTI I PROCESSORI      
// | DEVONO STAMPARE LA SOMMA E I TEMPI PASSANDO COME ULTIMO PARAMETRO -1.                            
// ----------------------------------------------------------------------------------------------------
// | ISTRUZIONI PER L'USO                                                                             
// | IL PRIMO ARGOMENTO DA PASSARE QUANDO SI DEVE ESEGUIRE IL COMANDO mpiexec E' LA QUANTITA' DI      
// | DI NUMERI DA SOMMARE E DEVE ESSERE UN NUMERO INTERO. SE QUESTO PARAMETRO E' <= 20 ALLORA BISOGNA 
// | INSERIRE TANTI NUMERI REALI QUANTI SPECIFICATI NEL PRIMO PARAMETRO. CIOE' SE SI VOGLIONO SOMMARE 
// | 3 NUMERI, ALLORA BISOGNA INSERIRE TRE NUMERI REALI DOPO IL PRIMO PARAMETRO. SE IL PRIMO PARAMETRO
// | > 20, ALLORA I NUMERI VERRANNO GENERATO CASUALMENTE E NON BISOGNA INSERIRE ALCUN NUMERO. IL      
// | PARAMETRO SUCCESSIVO DA INSERIRE E' QUELLO DELLA STRATEGIA DA USARE E DEVE ESSERE >= 1 E <= 3.   
// | L'ULTIMO PARAMETRO DA INSERIRE E QUELLO CHE INDICA IL PROCESSORE CHE STAMPA L'OUTPUT. DEVE       
// | ESSERE UN INTERO => -1 E <= NUMERO_PROCESSORI. SE SI VUOLE CHE TUTTI I PROCESSORI STAMPINO IL    
// | RISULTATO OTTENUTO BISOGNA INSERIRE -1, ALTRIMENTI SOLO IL PROCESSORE CON ID UGUALE AL PARAMETRO 
// | SPECIFICATO STAMPERA' IL RISULTATO                                                               
// ----------------------------------------------------------------------------------------------------
// | ESEMPI D'USO                                                                                     
// | CON IL SEGUENTE INPUT IL PROGRAMMA ESEGUIRA' LA SOMMA DI 5 NUMERI, OVVERO: 1, 2, 3, 4, 5, CON LA
// | STRATEGIA 1, CON 4 PROCESSORI E CON L'OPZIONE CHE IL PROCESSORE 2 STAMPERA' IL RISULTATO 
// | INPUT INSERITO: 
// | #!/bin/bash
// | #PBS -q studenti
// | #PBS -l nodes=8:ppn=8
// | #PBS -N numbers_5_strategy_1_printer_2
// | #PBS -o numbers_5_strategy_1_printer_2.out
// | #PBS -e numbers_5_strategy_1_printer_2.err
// | sort -u $PBS_NODEFILE > hostlist
// | NCPU=`wc -l < hostlist`
// | echo ------------------------------------------------------
// | echo ' This job is allocated on '${NCPU}' cpu(s)'
// | echo 'Job is running on node(s): '
// | cat hostlist
// | PBS_O_WORKDIR=$PBS_O_HOME/project_sum_of_numbers
// | echo ------------------------------------------------------
// | echo PBS: qsub is running on $PBS_O_HOST
// | echo PBS: originating queue is $PBS_O_QUEUE
// | echo PBS: executing queue is $PBS_QUEUE
// | echo PBS: working directory is $PBS_O_WORKDIR
// | echo PBS: execution mode is $PBS_ENVIRONMENT
// | echo PBS: job identifier is $PBS_JOBID
// | echo PBS: job name is $PBS_JOBNAME
// | echo PBS: node file is $PBS_NODEFILE
// | echo PBS: current home directory is $PBS_O_HOME
// | echo PBS: PATH = $PBS_O_PATH
// | echo ------------------------------------------------------
// | echo "Eseguo/usr/lib64/openmpi/1.4-gcc/bin/mpicc -o $PBS_O_WORKDIR/project_sum_of_numbers $PBS_O_WORKDIR/main.c"
// | /usr/lib64/openmpi/1.4-gcc/bin/mpicc -o $PBS_O_WORKDIR/project_sum_of_numbers $PBS_O_WORKDIR/main.c
// | echo "Eseguo:/usr/lib64/openmpi/1.4-gcc/bin/-machinefile $PBS_NODEFILE -n 2 $PBS_O_WORKDIR/project_sum_of_numbers"
// | /usr/lib64/openmpi/1.4-gcc/bin/mpiexec -machinefile hostlist -n 4 $PBS_O_WORKDIR/project_sum_of_numbers 5 1 2 3 4 5 1 2
// | OUTPUT GENERATO:
// ------------------------------------------------------
// | This job is allocated on 8 cpu(s)
// | Job is running on node(s): 
// | wn273.scope.unina.it
// | wn274.scope.unina.it
// | wn275.scope.unina.it
// | wn276.scope.unina.it
// | wn277.scope.unina.it
// | wn278.scope.unina.it
// | wn279.scope.unina.it
// | wn280.scope.unina.it
// ------------------------------------------------------
// | PBS: qsub is running on ui-studenti.scope.unina.it
// | PBS: originating queue is studenti
// | PBS: executing queue is studenti
// | PBS: working directory is /homes/DMA/PDC/2022/PCCGPP98C/project_sum_of_numbers
// | PBS: execution mode is PBS_BATCH
// | PBS: job identifier is 4004997.torque02.scope.unina.it
// | PBS: job name is numbers_5_strategy_1_printer_2
// | PBS: node file is /var/spool/pbs/aux//4004997.torque02.scope.unina.it
// | PBS: current home directory is /homes/DMA/PDC/2022/PCCGPP98C
// | PBS: PATH = /usr/lib64/openmpi/1.2.7-gcc/bin:/usr/kerberos/bin:/opt/exp_soft/unina.it/intel/composer_xe_2013_sp1.3.174/bin/intel64:/opt/exp_soft/unina.it/intel/composer_xe_2013_sp1.3.174/mpirt/bin/intel64:/opt/exp_soft/unina.it/intel/composer_xe_2013_sp1.3.174/bin/intel64:/opt/exp_soft/unina.it/intel/composer_xe_2013_sp1.3.174/bin/intel64_mic:/opt/exp_soft/unina.it/intel/composer_xe_2013_sp1.3.174/debugger/gui/intel64:/opt/d-cache/srm/bin:/opt/d-cache/dcap/bin:/opt/edg/bin:/opt/glite/bin:/opt/globus/bin:/opt/lcg/bin:/usr/local/bin:/bin:/usr/bin:/opt/exp_soft/HADOOP/hadoop-1.0.3/bin:/opt/exp_soft/unina.it/intel/composerxe/bin/intel64/:/opt/exp_soft/unina.it/MPJExpress/mpj-v0_38/bin:/homes/DMA/PDC/2022/PCCGPP98C/bin
// ------------------------------------------------------
// | Eseguo/usr/lib64/openmpi/1.4-gcc/bin/mpicc -o /homes/DMA/PDC/2022/PCCGPP98C/project_sum_of_numbers/project_sum_of_numbers /homes/DMA/PDC/2022/PCCGPP98C/project_sum_of_numbers/main.c
// | Eseguo:/usr/lib64/openmpi/1.4-gcc/bin/-machinefile /var/spool/pbs/aux//4004997.torque02.scope.unina.it -n 2 /homes/DMA/PDC/2022/PCCGPP98C/project_sum_of_numbers/project_sum_of_numbers
// | I'm the processor 2 - This is the sum 15.000000 - Total time elapsed: 0.010694 seconds
// ------------------------------------------------------
// | CON IL SEGUENTE INPUT IL PROGRAMMA ESEGUIRA' LA SOMMA DI 21 NUMERI CASUALE GENERATI 
// | AUTOMATICAMENTE CON LA STRATEGIA 2 CON 8 PROCESSORI E CON L'OPZIONE CHE TUTTI I PROCESSORI 
// | STAMPANO IL RISULTATO OVVERO INSERENDO -1 COME ULTIMO PARAMETRO
// | INPUT INSERITO:
// | #!/bin/bash
// | #PBS -q studenti
// | #PBS -l nodes=8:ppn=8
// | #PBS -N numbers_21_strategy_2_printer_minus_1
// | #PBS -o numbers_21_strategy_2_printer_minus_1.out
// | #PBS -e numbers_21_strategy_2_printer_minus_1.err
// | sort -u $PBS_NODEFILE > hostlist
// | NCPU=`wc -l < hostlist`
// | echo ------------------------------------------------------
// | echo ' This job is allocated on '${NCPU}' cpu(s)'
// | echo 'Job is running on node(s): '
// | cat hostlist
// | PBS_O_WORKDIR=$PBS_O_HOME/project_sum_of_numbers
// | echo ------------------------------------------------------
// | echo PBS: qsub is running on $PBS_O_HOST
// | echo PBS: originating queue is $PBS_O_QUEUE
// | echo PBS: executing queue is $PBS_QUEUE
// | echo PBS: working directory is $PBS_O_WORKDIR
// | echo PBS: execution mode is $PBS_ENVIRONMENT
// | echo PBS: job identifier is $PBS_JOBID
// | echo PBS: job name is $PBS_JOBNAME
// | echo PBS: node file is $PBS_NODEFILE
// | echo PBS: current home directory is $PBS_O_HOME
// | echo PBS: PATH = $PBS_O_PATH
// | echo ------------------------------------------------------
// | echo "Eseguo/usr/lib64/openmpi/1.4-gcc/bin/mpicc -o $PBS_O_WORKDIR/project_sum_of_numbers $PBS_O_WORKDIR/main.c"
// | /usr/lib64/openmpi/1.4-gcc/bin/mpicc -o $PBS_O_WORKDIR/project_sum_of_numbers $PBS_O_WORKDIR/main.c
// | echo "Eseguo:/usr/lib64/openmpi/1.4-gcc/bin/-machinefile $PBS_NODEFILE -n 2 $PBS_O_WORKDIR/project_sum_of_numbers"
// | /usr/lib64/openmpi/1.4-gcc/bin/mpiexec -machinefile hostlist -n 4 $PBS_O_WORKDIR/project_sum_of_numbers 21 2 -1
// ------------------------------------------------------
// | This job is allocated on 8 cpu(s)
// | Job is running on node(s): 
// | wn273.scope.unina.it
// | wn274.scope.unina.it
// | wn275.scope.unina.it
// | wn276.scope.unina.it
// | wn277.scope.unina.it
// | wn278.scope.unina.it
// | wn279.scope.unina.it
// | wn280.scope.unina.it
// ------------------------------------------------------
// | PBS: qsub is running on ui-studenti.scope.unina.it
// | PBS: originating queue is studenti
// | PBS: executing queue is studenti
// | PBS: working directory is /homes/DMA/PDC/2022/PCCGPP98C/project_sum_of_numbers
// | PBS: execution mode is PBS_BATCH
// | PBS: job identifier is 4004998.torque02.scope.unina.it
// | PBS: job name is numbers_21_strategy_2_printer_minus_1
// | PBS: node file is /var/spool/pbs/aux//4004998.torque02.scope.unina.it
// | PBS: current home directory is /homes/DMA/PDC/2022/PCCGPP98C
// | PBS: PATH = /usr/lib64/openmpi/1.2.7-gcc/bin:/usr/kerberos/bin:/opt/exp_soft/unina.it/intel/composer_xe_2013_sp1.3.174/bin/intel64:/opt/exp_soft/unina.it/intel/composer_xe_2013_sp1.3.174/mpirt/bin/intel64:/opt/exp_soft/unina.it/intel/composer_xe_2013_sp1.3.174/bin/intel64:/opt/exp_soft/unina.it/intel/composer_xe_2013_sp1.3.174/bin/intel64_mic:/opt/exp_soft/unina.it/intel/composer_xe_2013_sp1.3.174/debugger/gui/intel64:/opt/d-cache/srm/bin:/opt/d-cache/dcap/bin:/opt/edg/bin:/opt/glite/bin:/opt/globus/bin:/opt/lcg/bin:/usr/local/bin:/bin:/usr/bin:/opt/exp_soft/HADOOP/hadoop-1.0.3/bin:/opt/exp_soft/unina.it/intel/composerxe/bin/intel64/:/opt/exp_soft/unina.it/MPJExpress/mpj-v0_38/bin:/homes/DMA/PDC/2022/PCCGPP98C/bin
// ------------------------------------------------------
// | Eseguo/usr/lib64/openmpi/1.4-gcc/bin/mpicc -o /homes/DMA/PDC/2022/PCCGPP98C/project_sum_of_numbers/project_sum_of_numbers /homes/DMA/PDC/2022/PCCGPP98C/project_sum_of_numbers/main.c
// | Eseguo:/usr/lib64/openmpi/1.4-gcc/bin/-machinefile /var/spool/pbs/aux//4004998.torque02.scope.unina.it -n 2 /homes/DMA/PDC/2022/PCCGPP98C/project_sum_of_numbers/project_sum_of_numbers
// | I'm the processor 0 - This is the sum 358.887848 - Total time elapsed: 0.000057 seconds
// | I'm the processor 2 - This is the sum 358.887848 - Total time elapsed: 0.000057 seconds
// | I'm the processor 1 - This is the sum 358.887848 - Total time elapsed: 0.000057 seconds
// | I'm the processor 3 - This is the sum 358.887848 - Total time elapsed: 0.000057 seconds

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
	{
		float scale = rand() / (float) RAND_MAX;
		numbers_to_sum[i] = -100 + scale * (100 - (-100));
	}
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
	return number > 0 && ((number & (number-1)) == 0 || number == 1);
}

int get_tree_hight(int processors)
{
	return (int)(log10(processors)/log10(2));
}

void strategy1(float *sum, int id_processor, int processors, int id_printer_processor)
{
	if (id_processor == id_printer_processor)
	{
		int i;
		for(i = 0; i < processors; ++i)
		{
			if(i != id_printer_processor)
			{
				MPI_Status status;
				int id_comunication = 80 + i;
				float partial_sum = 0;
				MPI_Recv(&partial_sum, 1, MPI_FLOAT, i, id_comunication, MPI_COMM_WORLD, &status);
				*sum += partial_sum;
			}
		}
	}
	else
	{
		int id_comunication = id_processor + 80;
		MPI_Send(sum, 1, MPI_FLOAT, id_printer_processor, id_comunication, MPI_COMM_WORLD);
	}
}

int rotate_tree(int id_processor, int id_printer, int processors)
{
	return (processors - id_printer + id_processor) % processors;
}

void strategy2(float *sum, int id_processor, int processors, int id_printer)
{
	int i;
	for(i = 0; i < get_tree_hight(processors); ++i)
	{
		if(rotate_tree(id_processor, id_printer, processors) % (int)pow(2, i) == 0)
		{
			if(rotate_tree(id_processor, id_printer, processors) % (int)pow(2, i + 1) == 0)
			{
				int id_comunication = id_processor;
				float partial_sum = 0;
				MPI_Status status;
				int sender = (int)(id_processor + pow(2, i)) % processors;
				MPI_Recv(&partial_sum, 1, MPI_FLOAT, sender, id_comunication, MPI_COMM_WORLD, &status);
				*sum += partial_sum;				
			}
			else
			{	
				int destination = (((int)(id_processor - pow(2, i)) % processors) + processors) % processors;
				int id_comunication = destination;
				MPI_Send(sum, 1, MPI_FLOAT, destination, id_comunication, MPI_COMM_WORLD);
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

void sum_numbers(int strategy, float *sum, int size_numbers_to_sum, int id_processor, int processors, int id_printer)
{
	if(strategy == 1 || !is_a_power_of_2(processors))
	{
		if(id_processor == 0 && strategy != 1)
			puts("I'm using the strategy 1. The number of processors is not a power of 2");
		strategy1(sum, id_processor, processors, (id_printer == -1) ? 0 : id_printer);
	}
	else if(strategy == 2)
		strategy2(sum, id_processor, processors, (id_printer == -1) ? 0 : id_printer);
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

void send_sum_if_broadcast_true(int id_printer_processor, int id_processor, float *sum)
{
	if(id_printer_processor == -1)
		MPI_Bcast(sum, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
}

void send_total_time_elapsed_if_broadcast_true(int id_printer_processor, int id_processor, double *total_time_elapsed)
{
	if(id_printer_processor == -1)
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
	MPI_Bcast(&id_printer_processor, 1, MPI_INT, 0, MPI_COMM_WORLD);
	sum_numbers(strategy, &sum, size_numbers_to_sum, id_processor, processors, id_printer_processor);
	double t1 = MPI_Wtime();
	double time = t1 - t0;
	send_sum_if_broadcast_true(id_printer_processor, id_processor, &sum);
	double total_time_elapsed = 0;
	MPI_Reduce(&time, &total_time_elapsed, 1, MPI_DOUBLE, MPI_MAX, (id_printer_processor == -1) ? 0 : id_printer_processor, MPI_COMM_WORLD);
	send_total_time_elapsed_if_broadcast_true(id_printer_processor, id_processor, &total_time_elapsed);
	print_sum(id_processor, id_printer_processor, sum, total_time_elapsed);
	free_memory(id_processor, numbers_to_sum, local_numbers_to_sum);
	MPI_Finalize();
	return 0;
}
