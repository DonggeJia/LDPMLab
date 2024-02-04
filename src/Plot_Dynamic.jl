#using Plots
function Plot_Dynamic(displace, Boun, external)
inc = size(displace,2) # Number of Converged Steps to Plot
delta = 200;
a = findall(Boun.==3)  # DOFs Associated with Imposed Displacement
imp = [(Boun.==3); Bool.(zeros(length(nodes_of_steel_coor) * 6))]
steps = collect(1:delta:inc)
x = zeros(size(steps))
y = zeros(size(steps))
for i in eachindex(steps)
    step = steps[i]
    x[i] = displace[a[1],step].*sign(displace[a[1],end])
    y[i] = sum(external[imp,step]) .*sign(displace[a[1],end]) ./1000
end
return x, y
end
x, y = Plot_Dynamic(displace, Boun, external)
plot(scatter(x=x, y=y, mode="lines"))

# xlabel!("Displacement [mm]")
# ylabel!("Force [kN]")
