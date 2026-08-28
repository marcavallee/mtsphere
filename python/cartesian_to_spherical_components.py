import numpy as np

def cartesian_to_spherical_components(x, y, z, F):
    """
    Calcule les composantes sphériques (F_theta, F_phi) d'un champ vectoriel.
    Accepte des nombres isolés ou des tableaux NumPy.
    """
    Fx = F[:,:,0]
    Fy = F[:,:,1]
    Fz = F[:,:,2]

    # 1. Calcul des rayons et distances intermédiaires
    r_xy_sq = x ** 2 + y ** 2
    r_xy = np.sqrt(r_xy_sq)
    r = np.sqrt(r_xy_sq + z ** 2)

    # 2. Gestion de la singularité à l'origine (r = 0)
    # Évite les divisions par zéro en remplaçant temporairement les zéros
    is_origin = (r == 0)
    r_safe = np.where(is_origin, 1.0, r)
    r_xy_safe = np.where(r_xy == 0, 1.0, r_xy)

    # 3. Calcul des fonctions trigonométriques directes
    cos_theta = z / r_safe
    sin_theta = r_xy / r_safe

    cos_phi = x / r_xy_safe
    sin_phi = y / r_xy_safe

    # Cas particulier : sur l'axe Z (r_xy = 0), on fixe phi = 0 par convention
    cos_phi = np.where(r_xy == 0, 1.0, cos_phi)
    sin_phi = np.where(r_xy == 0, 0.0, sin_phi)

    # 4. Projection du champ vectoriel
    F_theta = Fx * cos_theta * cos_phi + Fy * cos_theta * sin_phi - Fz * sin_theta
    F_phi = -Fx * sin_phi + Fy * cos_phi

    # Si on est à l'origine, le champ sphérique n'est pas défini (on met 0.0)
    F_theta = np.where(is_origin, 0.0, F_theta)
    F_phi = np.where(is_origin, 0.0, F_phi)

    return F_theta, F_phi


# --- EXEMPLE D'UTILISATION ---
if __name__ == "__main__":
    # Exemple avec un point unique (X, Y, Z)
    x, y, z = 1.0, 1.0, 1.0
    # Exemple de champ : F(x,y,z) = (-y, x, 0) -> Un vortex pur
    Fx, Fy, Fz = -y, x, 0.0

    F_theta, F_phi = cartesian_to_spherical_components(x, y, z, Fx, Fy, Fz)

    print("--- Point unique ---")
    print(f"F_theta : {F_theta:.4f}")
    print(f"F_phi   : {F_phi:.4f}")  # Devrait être positif car le vortex tourne selon +phi
