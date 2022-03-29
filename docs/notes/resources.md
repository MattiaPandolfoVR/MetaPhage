---
sort: 6
---

# Required resources

Each process has been labelled as "big_res", "med_res", or "low_res",
and the current defaults assume a minimum 16 cores and 64 Gb of available
memory.

To run the pipeline in with less resources, you can add a custom configuration.
For example, to run with 8 cores and ~64 Gb of memory:

```text  
// Define the maximum resources
params {
  max_cpus = 8
  max_memory = 60.GB
  max_time = 72.h
}

// Configure the labels
process {
   withLabel:big_res {
        cpus = 8
        memory = 60.GB
        time = 24.h
    }
  withLabel: med_res {
        cpus = 8   
        memory = 32.GB
        time = 24.h
    }
  withLabel: low_res {
        cpus = 2   
        memory = 16.GB
        time = 24.h
    }
}
```
