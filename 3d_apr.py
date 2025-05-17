import os
import numpy as np
import csv
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

# Lists to hold combined data from all directories
temperatures_all = []
densities_all = []
vsa_all = []

# Iterate directories named by density
for dirname in sorted(os.listdir(".")):
    if os.path.isdir(dirname):
        try:
            density = float(dirname)  # Convert directory name to density float
        except:
            continue  # Skip if name is not a number

        for dirname1 in sorted(os.listdir(dirname)):
            if os.path.isdir(dirname+"/"+dirname1):
                try:
                    t = float(dirname1)
                except:
                    continue  # Skip if name is not a number
                temp_path = os.path.join(dirname+"/"+dirname1, "temperature.xvg")
                kinetic_path = os.path.join(dirname+"/"+dirname1, "kinetic.xvg")

                if os.path.isfile(temp_path) and os.path.isfile(kinetic_path):
                    try:
                        temp_data = np.loadtxt(temp_path, comments=('@', '#'))
                        kinetic_data = np.loadtxt(kinetic_path, comments=('@', '#'))

                        temp_data=temp_data[len(temp_data)//2:]  # Take the second half
                        kinetic_data=kinetic_data[len(kinetic_data)//2:]  # Take the second half

                        temp_values = temp_data[:, 1]
                        kinetic_values = kinetic_data[:, 1]
                        vsa=2*kinetic_values/0.0399478

                        # Append data points
                        temperatures_all.extend(temp_values)
                        densities_all.extend([density] * len(temp_values))
                        vsa_all.extend(vsa)

                        print(f"Loaded data from density={density} with {len(temp_values)} points")

                    except Exception as e:
                        print(f"Error reading files in {dirname}: {e}")

# Convert to numpy arrays
temperatures_all = np.array(temperatures_all)
densities_all = np.array(densities_all)
vsa_all = np.array(vsa_all)

# Prepare features for polynomial regression
X = np.vstack((temperatures_all, densities_all)).T
poly = PolynomialFeatures(degree=2, include_bias=False)
X_poly = poly.fit_transform(X)

# Fit model
model = LinearRegression()
model.fit(X_poly, vsa_all)

# Extract coefficients
coeff_labels = poly.get_feature_names_out(["temperature", "density"])
coeff_values = np.append(model.coef_, model.intercept_)

# Save coefficients to CSV
with open("vsa_vs_temp_density_poly2_coefficients.csv", "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Feature", "Coefficient"])
    for label, value in zip(np.append(coeff_labels, "intercept"), coeff_values):
        writer.writerow([label, value])

# ------------------ 3D Plotting --------------------
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Grid for surface plot
temp_range = np.linspace(temperatures_all.min(), temperatures_all.max(), 50)
density_range = np.linspace(densities_all.min(), densities_all.max(), 50)
temp_grid, density_grid = np.meshgrid(temp_range, density_range)

# Flatten grid to use in prediction
grid_points = np.vstack([temp_grid.ravel(), density_grid.ravel()]).T
grid_poly = poly.transform(grid_points)
vsa_pred = model.predict(grid_poly).reshape(temp_grid.shape)

# Plot regression surface
ax.plot_surface(temp_grid, density_grid, vsa_pred, color='orange', alpha=0.5)

# Optionally add original data points
ax.scatter(temperatures_all[::100], densities_all[::100], vsa_all[::100], color='blue', s=10, label='Data')

ax.set_xlabel('Temperature T, K')
ax.set_ylabel('Density rho, кг/м^3')
ax.set_zlabel('Average velocity square <v^2>, (м/с)^2')
ax.set_title('Polynomial Fit: Average velocity square vs Temperature and Density')

plt.legend()
plt.tight_layout()
plt.show()
