# BSPM: 3D-Kinetic Plasmasphere Model

The BSPM is a state-of-the-art 3D-Kinetic model for simulating the plasmasphere. This repository encompasses both the core model and a RESTful API interface for easy interaction.

## Quickstart

Before diving into the repository structure, here's a brief on how to get started:

1. **Initialize the API**:
   ```./start-api.sh```


2. **Manage the BSPM API Service on Ubuntu**:
- Enable the service:
  ```
  sudo systemctl enable start-bpim-api.service
  ```

- Start the service:
  ```
  sudo systemctl start start-bpim-api.service
  ```

- Restart the service (if needed):
  ```
  sudo systemctl restart start-bpim-api
  ```

- Check the service status:
  ```
  sudo systemctl status start-bpim-api
  ```

## Repository Structure

### `bspm/` - Core Model Codebase

- `BSPM_May2023.py`: This original script, introduces a new optional parameter `--executionid`. By default, this is set to the current date in the "YYYY-MM-DD" format. It influences the directory path of the output results.

- `execute.sh`: A utility shell script responsible for loading necessary libraries, executing the BSPM model, and returning execution messages for the RESTful API.

### `API/` - RESTful Interface

This directory contains the implementation of the RESTful API, designed using the FastAPI framework, which interfaces with the BSPM model for execution and result retrieval.

### Service Management

- `start-api.sh`: A bash script to boot up the API.

- `start-bspm-api.service`: A systemd service descriptor to manage the BSPM API as a service on Ubuntu systems.

## Additional Resources

For an in-depth understanding of the BSPM model or instructions on utilizing the API, please delve into the specific directories for comprehensive documentation.

---

**Note**: Should you encounter any issues or have further inquiries, kindly consult the extended documentation or reach out to the designated project maintainers.
