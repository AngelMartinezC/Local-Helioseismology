import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model


data = open("Data/solar_model_interior.txt","r")

data = [i for i in data]

r_R, cc, rho, P, gamma_1, T = [], [], [], [], [], []
for i in range(5,len(data)):
    a = data[i].split()
    r_R.append(float(a[0]))
    cc.append(float(a[1]))
    rho.append(float(a[2]))
    #P.append(float(a[3]))
    #gamma_1.append(float(a[4]))
    #T.append(float(a[5]))


def plot(x=0,y=0):
	f = plt.figure(figsize=(7,4))
	plt.plot(x,y,'g',label=r'$\rho$: density ')
	plt.ylabel('[$g/cm^3$]')
	plt.xlabel('Solar Radius $r/R$')
	plt.grid()
	plt.legend()
	plt.show()

def linear(x,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0):
	y = A+B*x+C*x**2+D*x**3+E*x**4+F*x**5+G*x**6+H*x**7+I*x**8+J*x**9+K*x**10
	return y

def deriv(x):
	y = B+2*C*x+3*D*x**2+4*E*x**3+5*F*x**4+6*G*x**5+7*H*x**6+8*I*x**7+\
	9*J*x**8+10*K*x**9
	return y

den = np.array(rho)*(100)**3/1000
r = np.array(r_R)*6.955e8
c = np.array(cc)/100

den = den[::-1]
r = r[::-1]
c = c[::-1]-min(c)



gmodel = Model(linear)
result = gmodel.fit(den,x=r,A=10,B=10,C=10,D=10,E=10,F=10,G=10,H=10,I=10,J=10,K=10)
Ad = result.best_values['A']
Bd = result.best_values['B']
Cd = result.best_values['C']
Dd = result.best_values['D']
Ed = result.best_values['E']
Fd = result.best_values['F']
Gd = result.best_values['G']
Hd = result.best_values['H']
Id = result.best_values['I']
Jd = result.best_values['J']
Kd = result.best_values['K']

gmodel = Model(linear)
result = gmodel.fit(c,x=r,A=1000,B=10,C=10,D=10,E=10,F=10,G=10,H=10,I=10,J=10,K=10)
Ac = result.best_values['A']
Bc = result.best_values['B']
Cc = result.best_values['C']
Dc = result.best_values['D']
Ec = result.best_values['E']
Fc = result.best_values['F']
Gc = result.best_values['G']
Hc = result.best_values['H']
Ic = result.best_values['I']
Jc = result.best_values['J']
Kc = result.best_values['K']

denfit = linear(x=r,A=Ad,B=Bd,C=Cd,D=Dd,E=Ed,F=Fd,G=Gd,H=Hd,I=Id,J=Jd,K=Kd)
cnew = linear(x=r,A=Ac,B=Bc,C=Cc,D=Dc,E=Ec,F=Fc,G=Gc,H=Hc,I=Ic,J=Jc,K=Kc)
cfit = cnew-min(cnew)
#print(result.fit_report())

#plt.plot(r,deriv(r))
#plt.show()




def L(l):
	y = np.sqrt(l*(l+1))
	return y

def lamb(r,l,c):
	y = L(l)*c/(2*np.pi*r)
	return y


XX = np.linspace(0,6955e8,100)

plt.plot(r[1::],lamb(r[1::],1,c=cfit[1::]))
plt.plot(r[1::],lamb(r[1::],5,c=cfit[1::]))
plt.plot(r[1::],lamb(r[1::],20,c=cfit[1::]))
plt.plot(r[1::],lamb(r[1::],50,c=cfit[1::]))
plt.plot(r[1::],lamb(r[1::],100,c=cfit[1::]))
plt.ylim(0,6e-3)
plt.show()








exit()
plt.plot(r,linear(x=r,A=A,B=B,C=C,D=D,E=E,F=F,G=G,H=H,I=I,J=J,K=K))
plt.plot(r,den)
plt.show()




 






















