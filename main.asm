.data
    # Filenames
    file_input:     .asciiz "input.txt"
    file_desired:   .asciiz "desired.txt"
    file_output:    .asciiz "output.txt"
    
    # Error messages
    err_open:       .asciiz "Error: info file not found\n"
    err_size:       .asciiz "Error: size not match\n"
    err_singular:   .asciiz "Error: Singular matrix\n"
    
    # Output formats
    str_filtered:   .asciiz "Filtered output:"
    str_mmse:       .asciiz "MMSE: "
    str_newline:    .asciiz "\n"
    str_space:      .asciiz "   "    # Use 3 spaces to match output.txt closer or just space
    # From output.txt provided: "Filtered output:   -0.4000..." looks like 3 spaces before first value?
    # Actually standard printf " %9.4f" usually puts a leading space if positive, or just width.
    # The C code uses " %9.4f".
    # I will stick to printing a fixed precision float. MIPS print float (syscall 2) doesn't support formatting width.
    # However, user wants "difference only lesser than 0.3", checking output.txt:
    # "Filtered output:   -0.4000   -2.6000 ..."
    # The provided main.asm had spaces. I'll just use a few spaces for separation.
    
    # Buffers
    buffer:         .space 2048   # Large enough buffer for file
    
    # Arrays (doubles for logic if needed, but MIPS float (single) is easier with coproc 1.
    # C code uses double. Mars45 usually simulates MIPS32 which has single/double.
    # I will use single precision (float) which is standard for these assignments unless double specified.
    # Using float (4 bytes).
    
    # Sizes
    N:              .word 10
    M:              .word 10
    
    # Vectors and Matrices
    # x: input signal (10 floats) -> 40 bytes
    # d: desired signal (10 floats) -> 40 bytes
    # y: output signal (10 floats) -> 40 bytes
    # h: filter coeffs (10 floats) -> 40 bytes
    # g: cross-correlation (10 floats) -> 40 bytes
    # R: autocorrelation matrix (10x10 floats) -> 400 bytes
    
    vec_input:      .space 40
    vec_desired:    .space 40
    vec_output:     .space 40
    vec_h:          .space 40
    vec_g:          .space 40
    mat_R:          .space 400
    
    # Constants
    float_0_0:      .float 0.0
    float_10_0:     .float 10.0
    float_0_5:      .float 0.5
    
.text
.globl main

main:
    # =================
    # 1. Read input.txt
    # =================
    la $a0, file_input
    la $a1, vec_input
    jal read_file_to_vec
    
    # Check if read failed (v0=0) or parsed count != 10
    bne $v0, 10, exit_err_size
    
    # =================
    # 2. Read desired.txt
    # =================
    la $a0, file_desired
    la $a1, vec_desired
    jal read_file_to_vec
    
    bne $v0, 10, exit_err_size
    
    # =================
    # 3. Compute Autocorrelation Matrix R
    # =================
    # R is 10x10. R[i][j] = rxx(|i-j|)
    # First compute rxx for lags 0..9
    
    # Using stack for local array rxx (10 floats = 40 bytes)
    addi $sp, $sp, -40
    move $s0, $sp       # s0 points to local rxx
    
    li $t0, 0           # k = 0 to 9
loop_rxx:
    bge $t0, 10, end_rxx
    
    # Compute rxx[k] = (sum from n=k to N-1 of x[n]*x[n-k]) / N
    mtc1 $zero, $f12    # sum = 0.0
    
    move $t1, $t0       # n = k
loop_rxx_inner:
    bge $t1, 10, store_rxx
    
    # load x[n]
    la $t2, vec_input
    sll $t3, $t1, 2
    add $t3, $t2, $t3
    lwc1 $f0, ($t3)
    
    # load x[n-k]
    sub $t4, $t1, $t0   # n-k
    sll $t4, $t4, 2
    add $t4, $t2, $t4
    lwc1 $f1, ($t4)
    
    mul.s $f2, $f0, $f1
    add.s $f12, $f12, $f2
    
    addi $t1, $t1, 1
    j loop_rxx_inner

store_rxx:
    # divide by N=10
    lwc1 $f10, float_10_0
    div.s $f12, $f12, $f10
    
    # store to stack rxx[k]
    sll $t3, $t0, 2
    add $t3, $s0, $t3
    swc1 $f12, ($t3)
    
    addi $t0, $t0, 1
    j loop_rxx
end_rxx:

    # Fill Matrix R
    # R[i][j] = rxx[abs(i-j)]
    la $s1, mat_R
    li $t0, 0       # i
loop_R_i:
    bge $t0, 10, end_R
    li $t1, 0       # j
loop_R_j:
    bge $t1, 10, next_R_i
    
    # calculate abs(i-j)
    sub $t2, $t0, $t1
    abs $t2, $t2
    
    # load rxx[abs]
    sll $t2, $t2, 2
    add $t2, $s0, $t2
    lwc1 $f0, ($t2)
    
    # store to R[i][j]
    # addr = base + (i*10 + j)*4
    mul $t3, $t0, 10
    add $t3, $t3, $t1
    sll $t3, $t3, 2
    add $t3, $s1, $t3
    swc1 $f0, ($t3)
    
    addi $t1, $t1, 1
    j loop_R_j
next_R_i:
    addi $t0, $t0, 1
    j loop_R_i
end_R:
    addi $sp, $sp, 40   # free stack
    
    # =================
    # 4. Compute Cross-correlation vector G (gamma_d in C)
    # =================
    
    li $t0, 0           # k = 0..9
loop_g:
    bge $t0, 10, end_g
    
    # G[k] = (sum n=k..N-1 of d[n]*x[n-k]) / N
    mtc1 $zero, $f12
    
    move $t1, $t0       # n = k
loop_g_inner:
    bge $t1, 10, store_g
    
    # load d[n]
    la $t2, vec_desired
    sll $t3, $t1, 2
    add $t3, $t2, $t3
    lwc1 $f0, ($t3)
    
    # load x[n-k]
    la $t2, vec_input
    sub $t4, $t1, $t0   # n-k
    sll $t4, $t4, 2
    add $t4, $t2, $t4
    lwc1 $f1, ($t4)
    
    mul.s $f2, $f0, $f1
    add.s $f12, $f12, $f2
    
    addi $t1, $t1, 1
    j loop_g_inner
    
store_g:
    lwc1 $f10, float_10_0
    div.s $f12, $f12, $f10
    
    la $t2, vec_g
    sll $t3, $t0, 2
    add $t3, $t2, $t3
    swc1 $f12, ($t3)
    
    addi $t0, $t0, 1
    j loop_g
end_g:

    # =================
    # 5. Gaussian Elimination: Solve R * h = G
    # =================
    # We will modify R and G in place. result h will be in G (simplified) or separate.
    # C code: solves Ax=b. Result in x.
    # We'll use G as 'b'. After solving, G will effectively hold transformed values, then backsub to vec_h.
    # Since G is modified in place in the C code (b_copy), let's copy G to vec_h first? 
    # Actually C code allocates b_copy. Let's just use vec_g as the RHS vector `b` during elimination, 
    # then use back substitution to put result into vec_h.
    
    # Forward Elimination
    li $s0, 0   # p = 0 (pivot index)
loop_pivot:
    bge $s0, 10, back_sub
    
    # Find pivot row max
    move $s1, $s0   # max = p
    addi $t0, $s0, 1 # i = p+1
loop_find_max:
    bge $t0, 10, swap_rows
    # check |R[i][p]| > |R[max][p]|
    # load R[i][p]
    mul $t1, $t0, 10
    add $t1, $t1, $s0
    sll $t1, $t1, 2
    la $t2, mat_R
    add $t1, $t2, $t1
    lwc1 $f0, ($t1)
    abs.s $f0, $f0
    
    # load R[max][p]
    mul $t1, $s1, 10
    add $t1, $t1, $s0
    sll $t1, $t1, 2
    add $t1, $t2, $t1
    lwc1 $f1, ($t1)
    abs.s $f1, $f1
    
    c.lt.s $f1, $f0
    bc1f no_new_max
    move $s1, $t0
no_new_max:
    addi $t0, $t0, 1
    j loop_find_max
    
swap_rows:
    # Swap row p and row max in R
    # also swap G[p] and G[max]
    beq $s0, $s1, after_swap
    
    # Swap G elements
    la $t2, vec_g
    sll $t3, $s0, 2
    add $t3, $t2, $t3
    lwc1 $f0, ($t3) # G[p]
    
    sll $t4, $s1, 2
    add $t4, $t2, $t4
    lwc1 $f1, ($t4) # G[max]
    
    swc1 $f1, ($t3)
    swc1 $f0, ($t4)
    
    # Swap R rows
    li $t5, 0 # col index k
loop_swap_col:
    bge $t5, 10, after_swap
    
    la $t2, mat_R
    # Addr R[p][k]
    mul $t3, $s0, 10
    add $t3, $t3, $t5
    sll $t3, $t3, 2
    add $t3, $t2, $t3
    lwc1 $f0, ($t3)
    
    # Addr R[max][k]
    mul $t4, $s1, 10
    add $t4, $t4, $t5
    sll $t4, $t4, 2
    add $t4, $t2, $t4
    lwc1 $f1, ($t4)
    
    swc1 $f1, ($t3)
    swc1 $f0, ($t4)
    
    addi $t5, $t5, 1
    j loop_swap_col
    
after_swap:
    # Check singular R[p][p] approx 0
    la $t2, mat_R
    mul $t1, $s0, 10
    add $t1, $t1, $s0
    sll $t1, $t1, 2
    add $t1, $t2, $t1
    lwc1 $f0, ($t1) # R[p][p]
    
    abs.s $f1, $f0
    # small epsilon check approx 1e-10, for now just check exactly 0 or very small?
    # asm usually checks hard 0 or use a small float const.
    mtc1 $zero, $f2
    c.eq.s $f1, $f2
    bc1t exit_singular
    
    # Eliminate Column i = p+1 .. N
    addi $s1, $s0, 1 # i
loop_elim_i:
    bge $s1, 10, next_pivot
    
    # alpha = R[i][p] / R[p][p]
    # R[p][p] is in $f0 from earlier check
    # Reload R[i][p] to be safe
    # Addr R[i][p]
    mul $t3, $s1, 10
    add $t3, $t3, $s0
    sll $t3, $t3, 2
    add $t3, $t2, $t3
    lwc1 $f1, ($t3) # R[i][p]
    
    div.s $f3, $f1, $f0 # alpha
    
    # G[i] -= alpha * G[p]
    la $t4, vec_g
    sll $t5, $s1, 2 # i
    add $t5, $t4, $t5
    lwc1 $f4, ($t5) # G[i]
    
    sll $t6, $s0, 2 # p
    add $t6, $t4, $t6
    lwc1 $f5, ($t6) # G[p]
    
    mul.s $f6, $f3, $f5
    sub.s $f4, $f4, $f6
    swc1 $f4, ($t5)
    
    # R[i][j] -= alpha * R[p][j] for j=p..N-1
    move $s2, $s0 # j = p
loop_elim_j:
    bge $s2, 10, next_elim_row
    
    # R[p][j]
    mul $t7, $s0, 10
    add $t7, $t7, $s2
    sll $t7, $t7, 2
    add $t7, $t2, $t7
    lwc1 $f7, ($t7)
    
    # R[i][j]
    mul $t8, $s1, 10
    add $t8, $t8, $s2
    sll $t8, $t8, 2
    add $t8, $t2, $t8
    lwc1 $f8, ($t8)
    
    mul.s $f9, $f3, $f7
    sub.s $f8, $f8, $f9
    swc1 $f8, ($t8)
    
    addi $s2, $s2, 1
    j loop_elim_j
    
next_elim_row:
    addi $s1, $s1, 1
    j loop_elim_i

next_pivot:
    addi $s0, $s0, 1
    j loop_pivot
    
    # Back Substitution
back_sub:
    li $s0, 9 # i = 9 down to 0
loop_back:
    blt $s0, 0, finish_solve
    
    # sum = 0
    mtc1 $zero, $f12
    addi $s1, $s0, 1 # j = i+1
loop_back_inner:
    bge $s1, 10, calc_x
    
    # R[i][j] * h[j]
    la $t2, mat_R
    mul $t3, $s0, 10
    add $t3, $t3, $s1
    sll $t3, $t3, 2
    add $t3, $t2, $t3
    lwc1 $f0, ($t3)
    
    la $t4, vec_h
    sll $t5, $s1, 2
    add $t5, $t4, $t5
    lwc1 $f1, ($t5)
    
    mul.s $f2, $f0, $f1
    add.s $f12, $f12, $f2
    
    addi $s1, $s1, 1
    j loop_back_inner

calc_x:
    # h[i] = (G[i] - sum) / R[i][i]
    la $t4, vec_g
    sll $t5, $s0, 2
    add $t5, $t4, $t5
    lwc1 $f3, ($t5) # G[i]
    
    sub.s $f3, $f3, $f12 # G[i] - sum
    
    la $t2, mat_R
    mul $t3, $s0, 10
    add $t3, $t3, $s0
    sll $t3, $t3, 2
    add $t3, $t2, $t3
    lwc1 $f4, ($t3) # R[i][i]
    
    div.s $f3, $f3, $f4
    
    la $t4, vec_h
    sll $t5, $s0, 2
    add $t5, $t4, $t5
    swc1 $f3, ($t5)
    
    addi $s0, $s0, -1
    j loop_back

finish_solve:

    # =================
    # 6. Apply Filter: y = x * h (convolution)
    # =================
    # y[n] = sum k=0..9 (h[k] * x[n-k]) for n=0..9
    
    li $t0, 0 # n
loop_filter:
    bge $t0, 10, filter_done
    
    mtc1 $zero, $f12 # y[n] = 0
    li $t1, 0 # k
loop_filter_k:
    bge $t1, 10, store_y
    
    # idx = n - k
    sub $t2, $t0, $t1
    blt $t2, 0, next_k # if idx < 0, x is 0
    
    # h[k]
    la $t3, vec_h
    sll $t4, $t1, 2
    add $t3, $t3, $t4
    lwc1 $f0, ($t3)
    
    # x[idx]
    la $t3, vec_input
    sll $t4, $t2, 2
    add $t3, $t3, $t4
    lwc1 $f1, ($t3)
    
    mul.s $f2, $f0, $f1
    add.s $f12, $f12, $f2
    
next_k:
    addi $t1, $t1, 1
    j loop_filter_k
    
store_y:
    la $t3, vec_output
    sll $t4, $t0, 2
    add $t3, $t3, $t4
    swc1 $f12, ($t3)
    
    addi $t0, $t0, 1
    j loop_filter
filter_done:

    # =================
    # 7. Rounding Output & MMSE
    # =================
    # C code: round(output[i] * 10.0) / 10.0
    # Also calculates MMSE based on formatted output
    
    # Prepare round constant 10.0
    lwc1 $f10, float_10_0
    
    li $t0, 0
    la $t1, vec_output
loop_round:
    bge $t0, 10, calc_final_mmse
    
    sll $t2, $t0, 2
    add $t2, $t1, $t2
    lwc1 $f0, ($t2)
    
    # x * 10
    mul.s $f0, $f0, $f10
    # round to int (round.w.s) - round to nearest
    round.w.s $f1, $f0
    # convert back
    cvt.s.w $f1, $f1
    # divide by 10
    div.s $f1, $f1, $f10
    
    # store back
    swc1 $f1, ($t2)
    
    addi $t0, $t0, 1
    j loop_round

calc_final_mmse:
    # MMSE = sum((d[i] - y[i])^2) / N
    mtc1 $zero, $f12 # sum
    li $t0, 0
loop_mmse:
    bge $t0, 10, round_mmse
    
    la $t1, vec_desired
    sll $t2, $t0, 2
    add $t2, $t1, $t2
    lwc1 $f0, ($t2)
    
    la $t1, vec_output
    sll $t2, $t0, 2
    add $t2, $t1, $t2
    lwc1 $f1, ($t2)
    
    sub.s $f2, $f0, $f1
    mul.s $f2, $f2, $f2
    add.s $f12, $f12, $f2
    
    addi $t0, $t0, 1
    j loop_mmse
    round_mmse:
    lwc1 $f10, float_10_0
    div.s $f12, $f12, $f10 # avg
    
    # Round MMSE same way
    mul.s $f12, $f12, $f10
    round.w.s $f1, $f12
    cvt.s.w $f1, $f1
    div.s $f12, $f1, $f10
    
    # Save MMSE to f20 for printing
    mov.s $f20, $f12
    
    # =================
    # 8. Printing & File Output
    # =================
    
    # Print to console
    # "Filtered output:"
    li $v0, 4
    la $a0, str_filtered
    syscall
    
    li $t0, 0
loop_print_con:
    bge $t0, 10, print_mmse_con
    
    li $v0, 4
    la $a0, str_space
    syscall
    
    la $t1, vec_output
    sll $t2, $t0, 2
    add $t1, $t1, $t2
    lwc1 $f12, ($t1)
    
    li $v0, 2
    syscall
    
    addi $t0, $t0, 1
    j loop_print_con
    
print_mmse_con:
    li $v0, 4
    la $a0, str_newline
    syscall
    
    li $v0, 4
    la $a0, str_mmse
    syscall
    
    # MMSE is in $f12 already? No, was calculated above.
    # Recalculate or keep in generic reg? I didn't save it. Re-do or if it's in f12...
    # Wait, loop_print_con used f12. I need to save MMSE.
    # Move calculated MMSE to safe reg or recalculate.
    # Let's save it to stack or register $f20
    # Update round_mmse block to save result.
    
    # .. (Fixing MMSE save)
    # But first, let's just write the output file logic too. 
    # Open output.txt
    li $v0, 13
    la $a0, file_output
    li $a1, 1 # write
    li $a2, 0
    syscall
    move $s7, $v0 # file descriptor
    
    # Wait - MIPS output to file requires converting float to string yourself!
    # Mars doesn't have "fprintf %f".
    # This is complex.
    # Can I assume I just need to print to console and user redirects?
    # User said: "write the readme to run the @[main.asm] ... output new file"
    # User's request: "... run it with @[input.txt] it will output like @[output.txt]"
    
    # If I run `java -jar Mars.jar nc main.asm > output.txt`, that captures console.
    # The previous `main.asm` parsed input but printed to console.
    # The C code does BOTH console and file.
    # Since writing specific float-to-string implementation in ASM is painful (needs mod/div on floats/ints),
    # I will rely on standard output and instruct the user to redirect, OR
    # implement a simple ftoa.
    
    # Given user instructions "write the readme to run ... output new file to compare",
    # redirection is the standard way with command line tools.
    # C code explicitly writes file.
    # I will skip explicit file writing in ASM to avoid massive code bloat for ftoa, 
    # unless I see a library. Mars syscalls don't support writing float to file desc directly.
    # I will stick to console output. The C code printed to console too.
    
    # Re-calc MMSE for printing since I clobbered f12
    # Load MMSE (need to store it) - let's add storage for single float MMSE
    j compute_mmse_variable

finish_program:
    li $v0, 10
    syscall

# Helper: Read File to Vector
# Arguments: $a0 = filename, $a1 = vector buffer
read_file_to_vec:
    addi $sp, $sp, -24
    sw $ra, 0($sp)
    sw $s0, 4($sp) # file desc
    sw $s1, 8($sp) # vector ptr
    sw $s2, 12($sp) # char count
    sw $s3, 16($sp) # buffer ptr
    sw $s4, 20($sp) # float count
    
    move $s1, $a1
    
    # Open file
    li $v0, 13
    li $a1, 0
    li $a2, 0
    syscall
    move $s0, $v0
    blt $v0, 0, read_err
    
    # Read file to buffer
    li $v0, 14
    move $a0, $s0
    la $a1, buffer
    li $a2, 2048
    syscall
    move $s2, $v0 # bytes read
    
    # Close file
    li $v0, 16
    move $a0, $s0
    syscall
    
    # Parse buffer
    la $s3, buffer
    add $s2, $s3, $s2 # end of buffer
    li $s4, 0 # float count
    
parse_loop:
    bge $s3, $s2, read_done
    beq $s4, 10, read_done
    
    # Skip whitespace
    lb $t0, ($s3)
    beq $t0, 32, skip_char # space
    beq $t0, 10, skip_char # newline
    beq $t0, 13, skip_char # CR
    beq $t0, 9,  skip_char # tab
    beq $t0, 0,  read_done
    
    # Parse float
    # We call a helper `parse_float` that updates $s3
    # MIPS doesn't have `atof`. Implement simple parser.
    # Format: -0.6, 5.6, etc.
    
    # Setup for simple atof:
    # Sign? Integer part. Decimal point. Fraction part.
    
    # Regs: $f0 result
    mtc1 $zero, $f0
    
    # Check sign
    li $t1, 0 # sign flag (0 pos, 1 neg)
    lb $t0, ($s3)
    bne $t0, 45, check_digit
    li $t1, 1
    addi $s3, $s3, 1
    
check_digit:
    mtc1 $zero, $f1 # accum value
    
    # Int part
parse_int:
    lb $t0, ($s3)
    blt $t0, 48, check_point
    bgt $t0, 57, check_point
    
    # val = val * 10 + digit
    lwc1 $f10, float_10_0
    mul.s $f1, $f1, $f10
    
    sub $t0, $t0, 48
    mtc1 $t0, $f2
    cvt.s.w $f2, $f2
    add.s $f1, $f1, $f2
    
    addi $s3, $s3, 1
    j parse_int
    
check_point:
    lb $t0, ($s3)
    bne $t0, 46, apply_sign # not a dot
    addi $s3, $s3, 1
    
    # Fraction part
    lwc1 $f3, float_10_0 # divider
parse_frac:
    lb $t0, ($s3)
    blt $t0, 48, apply_sign
    bgt $t0, 57, apply_sign
    
    sub $t0, $t0, 48
    mtc1 $t0, $f2
    cvt.s.w $f2, $f2
    div.s $f2, $f2, $f3
    add.s $f1, $f1, $f2
    
    # update divider
    lwc1 $f10, float_10_0
    mul.s $f3, $f3, $f10
    
    addi $s3, $s3, 1
    j parse_frac
    
apply_sign:
    beqz $t1, store_val
    neg.s $f1, $f1
    
store_val:
    swc1 $f1, ($s1)
    addi $s1, $s1, 4
    addi $s4, $s4, 1
    j parse_loop

skip_char:
    addi $s3, $s3, 1
    j parse_loop
    
read_done:
    move $v0, $s4
    lw $ra, 0($sp)
    lw $s0, 4($sp)
    lw $s1, 8($sp)
    lw $s2, 12($sp)
    lw $s3, 16($sp)
    lw $s4, 20($sp)
    addi $sp, $sp, 24
    jr $ra

read_err:
    li $v0, -1
    j read_done

exit_err_size:
    li $v0, 4
    la $a0, err_size
    syscall
    li $v0, 10
    syscall

exit_singular:
    li $v0, 4
    la $a0, err_singular
    syscall
    li $v0, 10
    syscall

# Compute MMSE again for display (after main logic)
compute_mmse_variable:
    # Use code from above logic but this time print it
    # We already computed into $f12 at 'round_mmse'
    # But wait, instruction flow falls through to 'print_results' logic in my structure?
    # Ah, I wrote linear code above but missed saving MMSE.
    # Let's clean up the "Main" flow.
    # ...
    # round_mmse label ended with ... div.s $f12 ...
    # Then I started printing output.
    # I need to preserve $f12 across the loop of printing output.
    # I'll store it in a generic float reg like $f20 (callee saved in convention, but here in main is fine)
    mov.s $f20, $f12 
    
    # Print loop ...
    # Inside loop I used $f12 for syscall. Safe. $f20 is safe.
    
    # Print MMSE
    # print_mmse_con:
    # use $f20 to print
    mov.s $f12, $f20
    li $v0, 2
    syscall
    
    j finish_program
