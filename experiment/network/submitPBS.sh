#!/bin/sh -f

###############################################################
#                                                             #
#    Steve's shell script for submitting a parallel MPICH     #
#    job to the PBS queue using the qsub command.             #
#                                                             #
###############################################################

#     Remarks: A line beginning with # is a comment.
#	       A line beginning with #PBS is a PBS directive.
#              PBS directives must come first; any directives
#                 after the first executable statement are ignored.
#

#NJOB=$1
#CMD=$2
echo $NJOB
echo $CMD
echo $NJOB"_stdout"
##########################
#                        #
#   The PBS directives   #
#                        #
##########################

#          Set the name of the job (up to 15 characters, 
#          no blank spaces, start with alphanumeric character)

#PBS -N $NJOB

#          Specify the number of nodes requested and the
#          number of processors per node. 

#PBS -l nodes=1:ppn=12

#          By default, the standard output and error streams are sent
#          to files in the current working directory with names:
#              job_name.osequence_number  <-  output stream
#              job_name.esequence_number  <-  error stream
#          where job_name is the name of the job and sequence_number 
#          is the job number assigned when the job is submitted.
#          Use the directives below to change the files to which the
#          standard output and error streams are sent.

# This will save the name of the .sh script #PBS -o $JNAME_stdout
#PBS -o $NJOB"_stdout"
#PBS -e $NJOB"_stderr"

#          The directive below directs that the standard output and
#          error streams are to be merged, intermixed, as standard
#          output. 

#PBS -j oe

#          Specify the maximum cpu and wall clock time. The wall
#          clock time should take possible queue waiting time into
#          account.  Format:   hhhh:mm:ss   hours:minutes:seconds
#          Be sure to specify a reasonable value here.
#          If the job does not finish by the time reached,
#          the job is terminated.

#PBS -l     cput=6:00:00
#PBS -l walltime=6:00:00

#          Specify the queue.  The CMU cluster currently has three queues:
#          "blue" and "green" and "red". 

#PBS -q batch

#          Specify the maximum amount of physical memory required per process.
#          kb for kilobytes, mb for megabytes, gb for gigabytes.
#          Take some care in setting this value.  Setting it too large
#          can result in your job waiting in the queue for sufficient
#          resources to become available.

#	#PBS -l pmem=512mb

#          PBS can send informative email messages to you about the
#          status of your job.  Specify a string which consists of
#          either the single character "n" (no mail), or one or more
#          of the characters "a" (send mail when job is aborted),
#          "b" (send mail when job begins), and "e" (send mail when
#          job terminates).  The default is "a" if not specified.
#          You should also specify the email address to which the
#          message should be send via the -M option.

#PBS -m abe
#PBS -m ae

#PBS -M stamatiad.st@gmail.com

#          Declare the time after which the job is eligible for execution.
#          If you wish the job to be immediately eligible for execution,
#          comment out this directive.  If you wish to run at some time in 
#          future, the date-time argument format is
#                      [DD]hhmm
#          If the day DD is not specified, it will default to today if the
#          time hhmm is in the future, otherwise, it defaults to tomorrow.
#          If the day DD is specified as in the future, it defaults to the
#          current month, otherwise, it defaults to next month.

# #PBS -a 2215  commented out

#          Specify the priority for the job.  The priority argument must be
#          an integer between -1024 and +1023 inclusive.  The default is 0.

#  #PBS -p 0

#          Define the interval at which the job will be checkpointed,
#          if checkpointing is desired, in terms of an integer number
#          of minutes of CPU time.

#  #PBS -c c=2

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

NCPU=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`
PBS_JOBNAME=$NJOB
echo ------------------------------------------------------
echo ' This job is allocated on '${NCPU}' cpu(s)'
echo 'Job is running on node(s): '
cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: number of nodes is $NNODES
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------


cd $PBS_O_WORKDIR

# run the program

$CMD
#/usr/lib64/mpich2/bin/mpirun ../../mechanism_simple/x86_64/special -mpi -nobanner -c "PARALLEL=1" -c "SIMPLIFIED=1" -c "CLUSTER_ID=0" -c "EXPERIMENT=0" -c "EXPID=1" -c "AA=12" -c "VCLAMP=0.000000" -c "ISBINARY=1" final.hoc



