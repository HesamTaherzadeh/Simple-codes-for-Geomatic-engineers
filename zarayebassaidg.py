from sympy.physics.vector import *
import sympy as sym

a=sym.symbols('a')
b=sym.symbols('b')
var = []

var.extend(input("please enter your variables uv,tlr,etc: "))
print(var)
if var==['u', 'v']:
    u,v=sym.symbols(var[0]),sym.symbols(var[1])
elif var==['t','l','r']:
    t, l= sym.symbols(var[0]), sym.symbols(var[1])
    r=sym.symbols('r')
else:
    x, y = sym.symbols(var[0]), sym.symbols(var[1])


Eq=[]
for i in range(3):
    Eq.append(eval(input("enter the equations in order :")))

varp=[]


N = sym.physics.vector.ReferenceFrame('N')
if var==['u', 'v']:
     varp1=sym.diff(Eq[0],u)*N.x+sym.diff(Eq[1],u)*N.y+sym.diff(Eq[2],u)*N.z
     varp2=sym.diff(Eq[0],v)*N.x+sym.diff(Eq[1],v)*N.y+sym.diff(Eq[2],v)*N.z
elif var==['t','l','r']:
    varp1=(sym.diff(Eq[0],t)*N.x)+(sym.diff(Eq[1],t)*N.y)+(sym.diff(Eq[2],t)*N.z)
    varp2=(sym.diff(Eq[0], l) * N.x) + (sym.diff(Eq[1], l) * N.y)+ (sym.diff(Eq[2], l) * N.z)
else:
    varp1 = (sym.diff(Eq[0], x) * N.x) + (sym.diff(Eq[1], x) * N.y) + (sym.diff(Eq[2], x) * N.z)
    varp2 = (sym.diff(Eq[0], y) * N.x) + (sym.diff(Eq[1], y) * N.y) + (sym.diff(Eq[2], y) * N.z)



E,F,G=sym.physics.vector.dot(varp1,varp1),sym.physics.vector.dot(varp1,varp2),sym.physics.vector.dot(varp2,varp2)

print(f"the first derivative is :{varp1}")
print(f"the second derivative is :{varp2}")




print("\t\t the first coefs")
print(f"E is ={sym.simplify(E)}")
print("------------------------------------------------\n")
print(f"F is ={sym.simplify(F)}")
print("------------------------------------------------\n")
print(f"G is ={sym.simplify(G)}")
print("------------------------------------------------\n")


if input("Do you want the angle?   ")=='y':
    cangle=F/sym.sqrt(E*G)
    print("\ncos(a)=F/(E*G)**1/2\n")
    sym.pprint(sym.acos(cangle))

print("------------------------------------------------\n")
print("EG-F**2=\n")
sym.pprint(sym.sqrt(E*G-F**2))
print("------------------------------------------------\n")
print("\t\t the second coefs")
cross=sym.physics.vector.cross(varp1,varp2)
crossm=cross.magnitude()
n=cross/crossm

if var==['u', 'v']:
    varp11 = varp1.diff(u,N)
    varp12 = varp1.diff(v,N)
    varp22 = varp2.diff(v,N)
elif var==['t','l','r']:
    varp11=varp1.diff(t,N)
    varp12=varp1.diff(l,N)
    varp22=varp2.diff(l,N)
else:
    varp1 = (sym.diff(Eq[0], x) * N.x) + (sym.diff(Eq[1], x) * N.y) + (sym.diff(Eq[2], x) * N.z)
    varp2 = (sym.diff(Eq[0], y) * N.x) + (sym.diff(Eq[1], y) * N.y) + (sym.diff(Eq[2], y) * N.z)

print("------------------------------------------------\n")

L,M,N=sym.physics.vector.dot(varp11,n),sym.physics.vector.dot(varp12,n),sym.physics.vector.dot(varp22,n)

print(f"the cross product is:\t\n{cross}")
print("\n\n\n")
print(f"the cross product mag:\t\n{crossm}")
print("\n\n\n\n\n")
print("n is :")
sym.pprint(n)
print("\n\n\n\n\n")
print("L is =\n")
sym.pprint(sym.simplify(L))
print("M is =\n")
sym.pprint(sym.simplify(M))
print("N is =\n")
sym.pprint(sym.simplify(N))
print("------------------------------------------------\n")
#main curvature finders
print(f"**reminder**:\n")
print(f"E={E}\nG={G}\nF={F}\nL={L}\nM={M}\nN={N}\n")
print("**for main curvature**")
print("main curvature",end="=\n")
print("\nG*L+E*N-2*F*M",end="=\n")
num=G*L+E*N-2*F*M
sym.pprint(sym.simplify(num))
print("\n(((G*L+E*N-2*F*M)**2)-4*(E*G-F**2)*(L*N-M**2))**1/2",end="=\n")
num1=(((G*L+E*N-2*F*M)**2)-(4*(E*G-F**2)*(L*N-M**2)))**1/2
sym.pprint(sym.simplify(num1))
print("\n2(EG-F**2)",end="=\n")
denum=(2*(E*G-F**2))
sym.pprint(sym.simplify(denum))
print("\n")
denum,num,num1=sym.simplify(denum),sym.simplify(num),sym.simplify(num1)
print("the main curvature is ",end="=\n")
sym.pprint(sym.simplify((num-num1)/denum))
print("\n")
sym.pprint(sym.simplify((num+num1)/denum))
print("------------------------------------------------\n")
print("\n**for gaussian curvature**")
print("",end="=\n")
print("\nL*N-M**2",end="=\n")
numg=sym.simplify(L*N-M**2)
sym.pprint(sym.simplify(numg))
print("\n")
print("\n(EG-F**2)",end="=\n")
denumg=denum/2
sym.pprint(sym.simplify(denumg))
print("\n")
print("Gaussian",end="=\n")
sym.pprint(sym.simplify(numg/denumg))