import sympy as sp 

x, y = sp.symbols('x y')

#input function
func = input("Enter the function f(x,y): ")
fn = sp.sympify(func)

x0 = float(input("Enter the value of x0: "))
y0 = float(input("Enter the value of y0: "))
x_tar = float(input("Enter the target value of x: "))
tol = float(input("Enter the value of error tolerance: "))

derivatives = [fn]
taylor_sum = y0
order = 1

while True:
    prev = derivatives[-1]
    next_der = sp.diff(prev, x) + sp.diff(prev, y) * fn 
    derivatives.append(next_der)

    deri_val = next_der.subs({x:x0, y:y0})
    if deri_val == 0:
        order += 1
        continue

    k = order + 1
    rhs = (tol * sp.factorial(k)) / abs(deri_val)
    x_chk = x0 + rhs**( 1 / k )

    print(f"Order {k}: x_chk = {x_chk}, x_tar = {x_tar}")
    
    if (x_chk >= x_tar):
        x_curr = x_tar - x0
        result = y0
        for i in range(1, k+1):
            term = ((x_curr) ** i) / sp.factorial(i) * derivatives[i-1].subs({x:x0, y:y0})
            result += term
        print(f"\nApproximate value of y({x_tar}) is {result.evalf()} using Taylor Series of order {k}")
        break

    order += 1
