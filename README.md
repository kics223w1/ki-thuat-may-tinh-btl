# Wiener Filter Implementation

This project implements a Wiener filter for signal processing in both C and MIPS assembly.

## Files Structure

- `main.c` - C implementation of the Wiener filter
- `main.s` - MIPS assembly implementation of the Wiener filter
- `Mars45.jar` - MARS MIPS simulator (required to run the assembly version)
- `desired.txt` - Desired signal file (10 floating-point numbers)
- `input.txt` - Input signal file (10 floating-point numbers)
- `output.txt` - Generated output file
- `desired/` - Directory containing sample desired signal files
- `input/` - Directory containing sample input signal files
- `output/` - Directory containing expected output files for comparison

## Prerequisites

### For C Version
- GCC compiler

### For MIPS Assembly Version
- Java Runtime Environment (JRE) - required to run MARS simulator
- MARS MIPS Simulator (`Mars45.jar` - included in the project)

## Running the C Version

```bash
# Compile
gcc -o wiener main.c -lm

# Setup input files (do this once)
cp desired/desired19-44-21_11-Nov-25_10_10.txt desired.txt
cp input/input19-44-21_11-Nov-25_10_10_1.txt input.txt

# Run
./wiener
```

## Comparing Results

To verify the assembly implementation produces correct results, compare the output with expected files:

```bash
# Run and compare
cp input/input1.txt input.txt
./wiener

# Check expected output
cat output/output1.txt
```


```bash
# Run and compare
cp input/input2.txt input.txt
./wiener

# Check expected output
cat output/output2.txt
```

```bash
# Run and compare
cp input/input3.txt input.txt
./wiener

# Check expected output
cat output/output3.txt
```

```bash
# Run and compare
cp input/input4.txt input.txt
./wiener

# Check expected output
cat output/output4.txt
```

```bash
# Run and compare
cp input/input5.txt input.txt
./wiener

# Check expected output
cat output/output5.txt
```
# ki-thuat-may-tinh-btl
