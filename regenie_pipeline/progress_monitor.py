import subprocess
import datetime
import time

def get_step2_num():
    '''get the nubmer of step 2 jobs array still queuing'''
    # run squeue command in shell and get the output
    # ensure full job names and job IDs is printed
    output = subprocess.check_output('squeue -u $USER -o "%.18i %.9P %.100j %.8u %.2t %.10M %.6D %R"', shell=True).decode('utf-8')
    # split the output into lines
    lines = output.split('\n')
    num_lines = len(lines)
    # get the job IDs and names
    job_ids = [line.split()[0] for line in lines[1:num_lines-1]]
    job_names = [line.split()[2] for line in lines[1:num_lines-1]]
    job_states = [line.split()[4] for line in lines[1:num_lines-1]]
    # find the names of step2 jobs by their names and get their ids at the same time
    step2_job_ids = []
    step2_job_names = []
    step2_job_states = []
    for job_id, job_name, job_state in zip(job_ids, job_names, job_states):
        if 'step2' in job_name:
            step2_job_ids.append(job_id)
            step2_job_names.append(job_name)      
            step2_job_states.append(job_state)      
    # some of the will be in job array format, e.g. step2[1-22]
    # get the job IDs of those job arrays, and get their names at the same time
    step2_array_ids = []
    step2_array_names = []
    for step2_job_id, step2_job_name in zip(step2_job_ids, step2_job_names):
        if '[' in step2_job_id:
            step2_array_ids.append(step2_job_id)
            step2_array_names.append(step2_job_name)
    
    # beak them into individual job IDs, for eaxmple step2[1-22] will be broken into step2_1, step2_2, ..., step2_22
    result_job_names = []
    for step2_array_id, step2_array_name in zip(step2_array_ids, step2_array_names):
        step2_array_range = step2_array_id.split('[')[1].split(']')[0].split('-')
        step2_array_start = int(step2_array_range[0])
        step2_array_end = int(step2_array_range[1])
        for i in range(step2_array_start, step2_array_end+1):
            result_job_names.append(step2_array_name + '_' + str(i))
    # now include the job names that are not in job array format
    result_job_names.extend([job_name for job_name in step2_job_names if job_name not in step2_array_names])

    # which step2_job_states os in running state
    step2_run_jobs = [job_name for job_name, job_state in zip(step2_job_names, step2_job_states) if job_state == 'R']
    print(step2_run_jobs, flush=True)
    return [len(result_job_names), len(step2_run_jobs)]

step2_num_jobs_lst = []
time_stamps = []
# run get_step2_num() function every hour
# if the number of step 2 jobs array still queuing is 0, then exit the loop
while True:
    print("hi", flush=True)
    # record the time
    time_stamps.append(datetime.datetime.now())
    # record the number of step 2 jobs array still queuing
    step2_num, step2_run_num = get_step2_num()
    step2_num_jobs_lst.append(step2_num)
    if len(step2_num_jobs_lst) == 0:
        break
    else:
        print('The time is {}'.format(time_stamps[-1]), flush=True)
        print('There are still {} step 2 jobs array still queuing'.format(step2_num_jobs_lst[-1]), flush=True)
        print('There are {} step 2 jobs running'.format(step2_run_num), flush=True)
        print('------------------------------------------', flush=True)
        # get the rate of change of step 2 jobs array still queuing
        if len(step2_num_jobs_lst) > 1:
            time_diff = (time_stamps[-1] - time_stamps[0]).total_seconds() / 3600
            step2_num_jobs_rate = (step2_num_jobs_lst[-1] - step2_num_jobs_lst[0]) / time_diff            
            # print the rate of change of step 2 jobs array still queuing
            print('The rate of change of step 2 jobs array still queuing is {} per hour'.format(step2_num_jobs_rate), flush=True)
        # sleep for 1 hour
        time.sleep(3600)