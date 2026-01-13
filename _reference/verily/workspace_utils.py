# workspace_utils.py (place in workspace root or create as package)
  import os
  import json
  import subprocess
  from typing import Dict

  def setup_aou_env(verbose: bool = False) -> Dict[str, str]:
      """
      Set All of Us workspace environment variables using wb CLI.
      
      Call this once at the start of any notebook:
          from workspace_utils import setup_aou_env
          setup_aou_env()
      
      Returns:
          dict: Environment variables that were set
      """
      # Extract workspace info
      workspace = json.loads(
          subprocess.run(
              ["wb", "workspace", "describe", "--format=json"],
              capture_output=True, text=True, check=True
          ).stdout
      )

      # Extract resources
      resources = json.loads(
          subprocess.run(
              ["wb", "resource", "list", "--format=json"],
              capture_output=True, text=True, check=True
          ).stdout
      )

      # Set project
      os.environ["GOOGLE_CLOUD_PROJECT"] = workspace["googleProjectId"]

      # Set buckets and CDR
      os.environ["WORKSPACE_CDR"] = ""

      for r in resources:
          if r["resourceType"] == "GCS_BUCKET":
              if "temporary-workspace-bucket" in r["id"]:
                  os.environ["WORKSPACE_TEMP_BUCKET"] = f"gs://{r['bucketName']}"
              elif "workspace-bucket" in r["id"]:
                  os.environ["WORKSPACE_BUCKET"] = f"gs://{r['bucketName']}"

          elif r["resourceType"] in ["BQ_DATASET", "BIGQUERY_DATASET"]:
              if os.environ.get("WORKSPACE_CDR") == "":
                  os.environ["WORKSPACE_CDR"] = f"{r['projectId']}.{r['datasetId']}"

      env_vars = {
          "GOOGLE_CLOUD_PROJECT": os.environ.get("GOOGLE_CLOUD_PROJECT"),
          "WORKSPACE_BUCKET": os.environ.get("WORKSPACE_BUCKET"),
          "WORKSPACE_TEMP_BUCKET": os.environ.get("WORKSPACE_TEMP_BUCKET"),
          "WORKSPACE_CDR": os.environ.get("WORKSPACE_CDR")
      }

      if verbose:
          print("✅ Workspace environment variables set:")
          for key, val in env_vars.items():
              print(f"  {key} = {val}")

      return env_vars


  # Convenience function - even simpler
  def get_workspace_cdr() -> str:
      """Get CDR dataset, setting env vars if needed"""
      if "WORKSPACE_CDR" not in os.environ:
          setup_aou_env()
      return os.environ["WORKSPACE_CDR"]


  def get_workspace_bucket() -> str:
      """Get workspace bucket, setting env vars if needed"""
      if "WORKSPACE_BUCKET" not in os.environ:
          setup_aou_env()
      return os.environ["WORKSPACE_BUCKET"]

  Usage in any notebook:
  # Option 1: Explicit setup
  from workspace_utils import setup_aou_env
  setup_aou_env(verbose=True)

  # Now use variables
  CDR = os.environ['WORKSPACE_CDR']
  bucket = os.environ['WORKSPACE_BUCKET']

  # Option 2: Even simpler - auto-setup on first use
  from workspace_utils import get_workspace_cdr, get_workspace_bucket
  CDR = get_workspace_cdr()
  bucket = get_workspace_bucket()

  Even Better: IPython Startup Script

  For truly zero-setup, create an IPython startup script that runs automatically:

  # ~/.ipython/profile_default/startup/00-workspace.py
  import os
  import json
  import subprocess

  try:
      # Only run if wb CLI is available and not already set
      if "WORKSPACE_CDR" not in os.environ:
          workspace = json.loads(
              subprocess.run(
                  ["wb", "workspace", "describe", "--format=json"],
                  capture_output=True, text=True, check=True, timeout=5
              ).stdout
          )

          resources = json.loads(
              subprocess.run(
                  ["wb", "resource", "list", "--format=json"],
                  capture_output=True, text=True, check=True, timeout=5
              ).stdout
          )

          os.environ["GOOGLE_CLOUD_PROJECT"] = workspace["googleProjectId"]
          os.environ["WORKSPACE_CDR"] = ""

          for r in resources:
              if r["resourceType"] == "GCS_BUCKET":
                  if "temporary-workspace-bucket" in r["id"]:
                      os.environ["WORKSPACE_TEMP_BUCKET"] = f"gs://{r['bucketName']}"
                  elif "workspace-bucket" in r["id"]:
                      os.environ["WORKSPACE_BUCKET"] = f"gs://{r['bucketName']}"
              elif r["resourceType"] in ["BQ_DATASET", "BIGQUERY_DATASET"]:
                  if os.environ.get("WORKSPACE_CDR") == "":
                      os.environ["WORKSPACE_CDR"] = f"{r['projectId']}.{r['datasetId']}"

          print("✅ All of Us workspace variables loaded")

  except Exception as e:
      # Silently fail if not in AoU workspace
      pass
