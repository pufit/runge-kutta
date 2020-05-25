## Почему тема устойчивости важна? Зачем нужны неявные методы?

#### Системы уравнений  
Жесткая:
$$\begin{cases} y_1' = -2  y_1 + y_2 + 2  sin{t} \\
y_2' = 998 y_1 - 999  y_2 + 999  (\cos{t} - \sin{t})
\end{cases}$$

Нежесткая: 
$$\begin{cases}
y_1' = -2  y_1 + y_2 + 2  \sin{t}\\
y_2' = y_1 - 2  y_2 + 2  (\cos{t} - \sin{t}))
\end{cases}$$
Обе системы имеют одинаковое решение:
$$
\begin{cases}
y_1 = 2  e ^ {-t} + \sin{t}\\
y_2 = 2  e ^ {-t} + \cos{t}
\end{cases}
$$

#### 

#### Получаем решения с фиксированной точностью


![png](output_9_0.png)

## A-stability
Рассмотрим уравнение $y'(x) = k \cdot y(x)$, его решением является $y(x) = e^{k \cdot x} + С$. Возьмем шаг $h > 0$, $x_{0} = 0$, $y_{0} = 1$ и решим уравнение числено, тогда для любого метода Рунге-Кутты $y_{n+1} = R(z) * y_{n}$, где $z = h \cdot k$. $R(z)$ называется функцией устойчивости. Если $|R(z)| < 1$, то $\lim\limits_{n \to \infty} y_n = 0$. Множество $\{z \in C | |R(z)| \le 1\}$ называется областью устойчивости. Если область устойчивости включает в себя отрицательную часть комплексной плоскости $C_{-} = \{z \in C | Re(z) \le 0\}$, то метод называется `A-устойчивым`. Область устойчивости любого явного метода является компактом, то есть явный метод не может быть A-устойчивым.    


Для явного и неявного метода Эйлера легко вычислить $ R(z) $ самостоятельно:  

- Явный метод:
$$
y_{n+1} = y_n + h \cdot  f(x_n, y_n)\\
y_{n+1} = y_n + h \cdot k \cdot y_n \\  
y_{n+1} = y_n \cdot (1 + h \cdot k)\\
y_{n+1} = y_n \cdot (1 + z)\\   
R(z) = 1 + z
$$
- Неявный метод:
$$
y_{n+1} = y_n + h \cdot f(x_n, y_{n+1})\\  
y_{n+1} = y_n + h \cdot k \cdot y_{n+1}\\
y_{n+1} \cdot (1 - h \cdot k) = y_n\\  
y_{n+1} = y_n \cdot \frac1{1 - z}\\
R(z) = \frac1{1 - z}  
$$

Есть две формулы для расчета $R(z)$, которые работают для обоих видов методов и одна только для явных методов:  
>- `RungeKutta.R_i1(self, z)`  
>- `RungeKutta.R_i2(self, z)`  
>- `RungeKutta.R_e(self, z, p)` где `p` - порядок точности, для `RK4` `p = 4`




![png](output_17_0.png)

## Примеры
Посмотрим, как работают численные методы Эйлера на примере уравнения $y'(x) = k \cdot y(x)$, $k = -15$, $x_0 = 0$, $y_0 = 1$ при разных величинах шага $h$.


![png](output_25_0.png)


![png](output_26_0.png)


## B-stability
Возьмем нелинейное уравнение $y' = f(x, y)$. Пусть y(x), z(x) - любые два решения. Метод является B-устойчивым, если из условия $(f(x, y(x)) - f(x, z(x)), y(x) - z(x)) \le 0$ следует, что $\forall h \ge 0 \;\; ||\hat{y_1} - y_1|| \le ||\hat{y_0} - y_0||$, где $\hat{y_1}$ и $y_1$ численные аппроксимации  после одного шага при начальных условиях $\hat{y_0}$ и $y_0$ соответственно. Если метод удовлетворяет двум критериям::
>- $b_{i} \ge 0, i = 1,\ldots,s$  
>- $M = (m_{ij}) = (b_{i} \cdot a_{ij} + b_{j} \cdot a_{ji} - b_{i} \cdot b_{j}), i = 1,\ldots,s, j = 1,\ldots,s$ неотрицательно определена  

то метод является `B-устойчивым`. Можно проверить с помощью метода `RungeKutta.is_B_stable(self)`.


```python
print('Is RK4 B-stable:', rk4.is_B_stable())
print('Is Radau IA B-stable:', ria.is_B_stable())
print('Is forward Euler method B-stable:', fe.is_B_stable())
print('Is backward Euler method B-stable:', be.is_B_stable())
```

    Is RK4 B-stable: no
    Is Radau IA B-stable: yes
    Is forward Euler method B-stable: no
    Is backward Euler method B-stable: yes


## L-stability
Если метод является `A-устойчивым` и $\lim\limits_{z \to \infty}R(z) = 0$, то метод является `L-устойчивым`. Если матрица `a` неявного метода не вырождена и удовлетворяет любому из двух критериев:  
>- $a_{sj} = b_{j}, j = 1,\ldots,s$
>- $a_{i1} = b_{i}, i = 1,\ldots,s$  

то метод является `L-устойчивым`. Можно проверить с помощью метода `RungeKutta.is_L_stable(self)`.


```python
print('Is RK4 L-stable:', rk4.is_L_stable())
print('Is Radau IA L-stable:', ria.is_L_stable())
print('Is forward Euler method L-stable:', fe.is_L_stable())
print('Is backward Euler method L-stable:', be.is_L_stable())
```

    Is RK4 L-stable: no
    Is Radau IA L-stable: maybe
    Is forward Euler method L-stable: no
    Is backward Euler method L-stable: yes


## Order Stars
Область определяется так: $\{z \in C \;\big| \;|R(z)| > |e^z|\}$. Есть много разных теорем, в качестве примера можно привести факт, для явных методов Рунге-Кутты число листов звезды всегда равно порядку точности метода + 1 и при $z\to0$ угловые размеры всех $2 \cdot (p + 1)$ секторов равны 


![png](output_32_0.png)


## Литература
Hairer, Ernst; Wanner, Gerhard (1991), Solving ordinary differential equations II: Stiff and differential-algebraic problems (стр. 40-59)  
Lambert, J.D (1991), Numerical Methods for Ordinary Differential Systems. The Initial Value Problem (стр. 232-237)  
Iserles, Arieh (1996), A First Course in the Numerical Analysis of Differential Equations (стр. 56-63)  
