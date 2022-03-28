---
sort: 2
---

# Metadata file

The required metadata file is a CSV file. The first column has to be `Sample`, and there
are no other requirements. Other columns can be used for plotting or diversity analysis.

## Example

The following metadata (subset from [Liang et al 2020](https://www.nature.com/articles/s41586-020-2192-1)):

Sample     | Infant_feeding | Infant_formula_type | Infant_delivery_type  | Infant_gender 
-----------|----------------|---------------------|----------------|---------------
Sample3218 | Formula        | cow-milk based      | C-Section      | Male          
Sample3221 | Formula        | cow-milk based      | C-Section      | Male          
Sample3248 | Formula        | cow-milk based      | C-Section      | Male          
Sample3247 | Formula        | soy-protein based   | C-Section      | Female        
Sample3084 | Formula        | cow-milk based      | Vaginal delivery | Female        
Sample2914 | Formula        | cow-milk based      | Vaginal delivery | Female        
Sample2969 | Formula        | cow-milk based      | Vaginal delivery | Male          

Can be recorded as CSV as follows:

```text
Sample,Infant_feeding,Infant_formula_type,Infant_delivery_type,Infant_gender
Sample3218,Formula,cow-milk based,C-Section,Male
Sample3221,Formula,cow-milk based,C-Section,Male
Sample3248,Formula,cow-milk based,C-Section,Male
Sample3247,Formula,soy-protein based,C-Section,Female
Sample3084,Formula,cow-milk based,Vaginal delivery,Female
Sample2914,Formula,cow-milk based,Vaginal delivery,Female
Sample2969,Formula,cow-milk based,Vaginal delivery,Male
```