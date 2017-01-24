def fn_a(var_a, var_b):
    return var_a + var_b

def fn_b(var_a, var_b):
    return var_a - var_b

fn_pointer = fn_a

print fn_pointer(1,2)