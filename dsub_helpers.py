"""
Dsub Helper Functions for All of Us Research
============================================

Utility functions for managing dsub jobs in the All of Us Researcher Workbench.
These functions simplify common dsub operations like job submission, monitoring, and cleanup.

Author: BW
Created: June 2025
Dependencies: os, subprocess, re

Usage:
    from utilities.dsub_helpers import dsub_script, check_dsub_status, cancel_job
    
    # Submit a job
    dsub_script(
        label='my-analysis',
        machine_type='n1-standard-4',
        envs={'COHORT': 'diabetes'},
        in_params={'INPUT_FILE': 'gs://bucket/data.csv'},
        out_params={'OUTPUT_FILE': 'gs://bucket/results.csv'}
    )
    
    # Check job status
    check_dsub_status(age='1d')
"""

import os
import subprocess
import re
from typing import Dict, Optional, Union

def validate_age_format(age: str) -> bool:
    """
    Validate age format for dsub dstat command
    
    Parameters:
    -----------
    age : str
        Age string to validate
        
    Returns:
    --------
    bool
        True if format is valid, False otherwise
        
    Examples:
    ---------
    >>> validate_age_format('3d')
    True
    >>> validate_age_format('12h')
    True
    >>> validate_age_format('invalid')
    False
    """
    # Pattern: one or more digits followed by exactly one valid unit
    pattern = r'^\d+[smhdw]$'
    return bool(re.match(pattern, age.lower()))

def dsub_script(
    label: str,
    machine_type: str,
    envs: Dict[str, str],
    in_params: Dict[str, str],
    out_params: Dict[str, str],
    boot_disk: int = 100,
    disk_size: int = 150,
    image: str = 'us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.14',
    script: str = 'pgen_subset_multiancestry.sh',
    preemptible: bool = True
) -> None:
    """
    Submit a dsub job with specified parameters
    
    Parameters:
    -----------
    label : str
        Job label/name for identification
    machine_type : str
        GCP machine type (e.g., 'n1-standard-4', 'c4-standard-8')
    envs : dict
        Environment variables to set in the job
    in_params : dict
        Input parameters mapping parameter names to file paths
    out_params : dict
        Output parameters mapping parameter names to output paths
    boot_disk : int, default 100
        Boot disk size in GB
    disk_size : int, default 150
        Additional disk size in GB
    image : str, default All of Us image
        Docker image to use for the job
    script : str, default 'pgen_subset_multiancestry.sh'
        Script file to execute
    preemptible : bool, default True
        Whether to use preemptible instances
        
    Examples:
    ---------
    >>> dsub_script(
    ...     label='bmi-analysis',
    ...     machine_type='n1-standard-4',
    ...     envs={'ANALYSIS_TYPE': 'gwas'},
    ...     in_params={'INPUT_DATA': 'gs://bucket/cohort.csv'},
    ...     out_params={'RESULTS': 'gs://bucket/results/'}
    ... )
    """
    # Get user information
    dsub_user_name = os.getenv("OWNER_EMAIL").split('@')[0]
    user_name = os.getenv("OWNER_EMAIL").split('@')[0].replace('.', '-')
    job_name = f'{label}'
    
    # Build dsub command
    dsub_cmd = 'dsub '
    dsub_cmd += '--provider google-cls-v2 '
    dsub_cmd += f'--machine-type "{machine_type}" '
    
    if preemptible:
        dsub_cmd += '--preemptible '
        
    # Use hyperdisk for c4 machines, SSD for others
    if 'c4' in machine_type:
        dsub_cmd += '--disk-type "hyperdisk-balanced" '
    else:
        dsub_cmd += '--disk-type "pd-ssd" '
        
    dsub_cmd += f'--boot-disk-size {boot_disk} '
    dsub_cmd += f'--disk-size {disk_size} '
    dsub_cmd += '--user-project "${GOOGLE_PROJECT}" '
    dsub_cmd += '--project "${GOOGLE_PROJECT}" '
    dsub_cmd += f'--image "{image}" '
    dsub_cmd += '--network "network" '
    dsub_cmd += '--subnetwork "subnetwork" '
    dsub_cmd += '--service-account "$(gcloud config get-value account)" '
    dsub_cmd += f'--user "{dsub_user_name}" '
    dsub_cmd += '--logging "${WORKSPACE_BUCKET}/dsub/logs/{job-name}/{user-id}/$(date +\'%Y%m%d\')/{job-id}-{task-id}-{task-attempt}.log" '
    dsub_cmd += ' "$@" '
    dsub_cmd += f'--name "{job_name}" '
    dsub_cmd += '--env GOOGLE_PROJECT="${GOOGLE_PROJECT}" '
    dsub_cmd += f'--script "{script}" '
    
    # Add environment variables
    for env_key, env_value in envs.items():
        dsub_cmd += f'--env {env_key}="{env_value}" '
        
    # Add input parameters
    for in_key, in_value in in_params.items():
        dsub_cmd += f'--input {in_key}="{in_value}" '
        
    # Add output parameters
    for out_key, out_value in out_params.items():
        dsub_cmd += f'--output {out_key}="{out_value}" '
    
    print(f"Submitting dsub job: {job_name}")
    print(f"Command: {dsub_cmd}")
    os.system(dsub_cmd)
    print("Job submitted successfully!")
    print()

def check_dsub_status(user: Optional[str] = None, full: bool = False, age: str = '1d') -> subprocess.CompletedProcess:
    """
    Check status of dsub jobs for the specified user
    
    Parameters:
    -----------
    user : str, optional
        Username to check jobs for. Defaults to current user from OWNER_EMAIL
    full : bool, default False
        Include full job details in output
    age : str, default '1d'
        Maximum age of jobs to display. Format: <integer><unit>
        Units: s (seconds), m (minutes), h (hours), d (days), w (weeks)
        Examples: '3d', '12h', '30m', '7w'
        
    Returns:
    --------
    subprocess.CompletedProcess
        Result of the dstat command
        
    Examples:
    ---------
    >>> check_dsub_status(age='3d', full=True)  # Last 3 days, full details
    >>> check_dsub_status()  # Default: last day, summary view
    """
    if user is None:
        # Get current user if not specified
        user = os.getenv("OWNER_EMAIL").split('@')[0]
    
    project = os.getenv("GOOGLE_PROJECT")
    
    # Validate age parameter
    if age is not None:
        if not validate_age_format(age):
            raise ValueError(
                f"Invalid age format: '{age}'. "
                "Expected format: <integer><unit> where unit is one of: s, m, h, d, w. "
                "Examples: '3d', '12h', '30m', '7w'"
            )
    
    # Build command
    cmd_parts = [
        "dstat",
        "--provider google-cls-v2",
        f"--user {user}",
        "--status '*'",
        f"--project {project}"
    ]
    
    if full:
        cmd_parts.append("--full")
    
    if age:
        cmd_parts.append(f"--age {age}")
    
    cmd = " ".join(cmd_parts)
    print(f"Running: {cmd}")
    return subprocess.run(cmd, shell=True, capture_output=False)

def job_details(user: Optional[str] = None, job: Optional[str] = None) -> subprocess.CompletedProcess:
    """
    List detailed information for jobs, including failed ones
    
    Parameters:
    -----------
    user : str, optional
        Username to check jobs for. Defaults to current user
    job : str, optional
        Specific job ID to check. If None, checks all jobs
        
    Returns:
    --------
    subprocess.CompletedProcess
        Result of the dstat command
        
    Examples:
    ---------
    >>> job_details()  # All jobs for current user
    >>> job_details(job='my-job-id')  # Specific job details
    """
    project = os.getenv("GOOGLE_PROJECT")
    
    if user is None:
        user = os.getenv("OWNER_EMAIL").split('@')[0]
        
    if job is None:
        job_filter = "'*' "
    else:
        job_filter = f'--jobs {job} '
    
    cmd = f"dstat --provider google-cls-v2 --project {project} --user {user} --status {job_filter}--full"
    print(f"Running: {cmd}")
    return subprocess.run(cmd, shell=True, capture_output=False)

def cancel_job(job_id: str) -> subprocess.CompletedProcess:
    """
    Cancel a specific dsub job
    
    Parameters:
    -----------
    job_id : str
        ID of the job to cancel
        
    Returns:
    --------
    subprocess.CompletedProcess
        Result of the ddel command
        
    Examples:
    ---------
    >>> cancel_job('my-job-12345')
    """
    project = os.getenv("GOOGLE_PROJECT")
    
    cmd = f"ddel --provider google-cls-v2 --project {project} --jobs {job_id}"
    print(f"Running: {cmd}")
    print(f"Canceling job: {job_id}")
    return subprocess.run(cmd, shell=True, capture_output=False)

def cancel_running_jobs(user: Optional[str] = None) -> subprocess.CompletedProcess:
    """
    Cancel all running/pending jobs for a user (use with caution!)
    
    Parameters:
    -----------
    user : str, optional
        Username to cancel jobs for. Defaults to current user
        
    Returns:
    --------
    subprocess.CompletedProcess
        Result of the ddel command
        
    Warning:
    --------
    This will cancel ALL running jobs for the specified user.
    Use with caution!
    
    Examples:
    ---------
    >>> cancel_running_jobs()  # Cancel all your running jobs
    """
    project = os.getenv("GOOGLE_PROJECT")
    
    if user is None:
        user = os.getenv("OWNER_EMAIL").split('@')[0]
    
    # Cancel running jobs
    cancel_cmd = f"ddel --provider google-cls-v2 --project {project} --users '{user}' --jobs '*'"
    print(f"WARNING: This will cancel ALL running jobs for user '{user}'")
    print(f"Running: {cancel_cmd}")
    
    return subprocess.run(cancel_cmd, shell=True, capture_output=False)

def get_job_logs(job_name: str, user: Optional[str] = None) -> str:
    """
    Get the log path for a specific job
    
    Parameters:
    -----------
    job_name : str
        Name of the job
    user : str, optional
        Username. Defaults to current user
        
    Returns:
    --------
    str
        Path to job logs in workspace bucket
    """
    if user is None:
        user = os.getenv("OWNER_EMAIL").split('@')[0]
    
    workspace_bucket = os.getenv("WORKSPACE_BUCKET")
    user_id = user.replace('.', '-')
    
    # Log path pattern matches what's used in dsub_script
    log_path = f"{workspace_bucket}/dsub/logs/{job_name}/{user_id}/"
    
    print(f"Job logs location: {log_path}")
    print("Use gsutil ls to see specific log files:")
    print(f"gsutil ls {log_path}*/")
    
    return log_path

# Convenience functions for common machine types
def submit_standard_job(label: str, script: str, envs: Dict[str, str] = None, 
                       in_params: Dict[str, str] = None, out_params: Dict[str, str] = None) -> None:
    """Submit job with standard n1-standard-4 machine"""
    dsub_script(
        label=label,
        machine_type='n1-standard-4',
        script=script,
        envs=envs or {},
        in_params=in_params or {},
        out_params=out_params or {}
    )

def submit_highmem_job(label: str, script: str, envs: Dict[str, str] = None,
                      in_params: Dict[str, str] = None, out_params: Dict[str, str] = None) -> None:
    """Submit job with high-memory n1-highmem-8 machine"""
    dsub_script(
        label=label,
        machine_type='n1-highmem-8',
        script=script,
        envs=envs or {},
        in_params=in_params or {},
        out_params=out_params or {},
        disk_size=200  # Higher disk for high-mem jobs
    )

def submit_compute_job(label: str, script: str, envs: Dict[str, str] = None,
                      in_params: Dict[str, str] = None, out_params: Dict[str, str] = None) -> None:
    """Submit job with compute-optimized c4-standard-8 machine"""
    dsub_script(
        label=label,
        machine_type='c4-standard-8',
        script=script,
        envs=envs or {},
        in_params=in_params or {},
        out_params=out_params or {},
        disk_size=300  # Higher disk for compute jobs
    )
