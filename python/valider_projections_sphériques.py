import numpy as np

def valider_projections_spheriques(X, Z, E_cartesian, E_theta, nom_champ="E"):
    """
    Effectue des tests de cohérence géométrique aux points cardinaux.

    Parameters:
    -----------
    X, Z : np.ndarray (2D)
        Grilles de coordonnées cartésiennes issues de np.meshgrid(x, z).
    E_cartesian : np.ndarray (3D)
        Champ cartésien original de forme (nx, nz, 3) contenant [Ex, Ey, Ez].
    E_theta : np.ndarray (2D)
        Votre matrice calculée de la composante theta de forme (nx, nz).
    nom_champ : str
        Nom pour l'affichage ("E" ou "H").
    """
    print(f"=== RAPPORT DE VALIDATION GÉOMÉTRIQUE : CHAMP {nom_champ} ===")

    # Séparation des composantes cartésiennes [nx, nz, 3]
    Fx = E_cartesian[:, :, 0]
    Fy = E_cartesian[:, :, 1]
    Fz = E_cartesian[:, :, 2]

    # Écart de tolérance numérique pour les comparaisons de float
    tol = 1e-5

    # -------------------------------------------------------------
    # TEST 1 : Pôle Nord (X = 0, Z maximum positif)
    # -------------------------------------------------------------
    # On cherche l'indice où X vaut 0 (au centre horizontal) et Z est maximal (haut)
    idx_center_x = np.argmin(np.abs(X[:, 0]))  # Ligne ou colonne centrale en X
    idx_max_z = np.argmax(Z[idx_center_x, :])  # Z max le long de cet axe

    # On extrait les valeurs numériques à cet indice
    val_F_theta_pole = E_theta[idx_center_x, idx_max_z]
    val_Fx_pole = Fx[idx_center_x, idx_max_z]

    diff_pole = np.abs(val_F_theta_pole - val_Fx_pole)
    status_pole = "✓ SUCCÈS" if diff_pole < tol else "❌ ÉCHEC"

    print(f"\n[Test 1] Pôle Nord géométrique (X=0, Z={Z[idx_center_x, idx_max_z]:.1f}m) :")
    print(f"  - {nom_champ}_theta calculé : {val_F_theta_pole:e}")
    print(f"  - {nom_champ}_x cartésien   : {val_Fx_pole:e}")
    print(f"  -> Statut : {status_pole} (Écart: {diff_pole:e})")
    print(f"  (Note: Au pôle Nord, thêta s'aligne géométriquement avec +X)")

    # -------------------------------------------------------------
    # TEST 2 : Équateur Droit (X maximum positif, Z = 0)
    # -------------------------------------------------------------
    # On cherche l'indice où Z vaut 0 et X est maximal (droite)
    idx_center_z = np.argmin(np.abs(Z[0, :]))  # Ligne ou colonne centrale en Z
    idx_max_x = np.argmax(X[:, idx_center_z])  # X max le long de cet axe

    val_F_theta_equateur = E_theta[idx_max_x, idx_center_z]
    val_Fz_equateur = Fz[idx_max_x, idx_center_z]

    # Théoriquement, F_theta à l'équateur droit vaut exactement -Fz.
    # Dans votre problème, puisque Fz = 0, F_theta doit être proche de 0.
    diff_equateur = np.abs(val_F_theta_equateur - (-val_Fz_equateur))
    status_equateur = "✓ SUCCÈS" if diff_equateur < tol else "❌ ÉCHEC"

    print(f"\n[Test 2] Équateur droit (X={X[idx_max_x, idx_center_z]:.1f}m, Z=0) :")
    print(f"  - {nom_champ}_theta calculé : {val_F_theta_equateur:e}")
    print(f"  - -{nom_champ}_z cartésien  : {-val_Fz_equateur:e}")
    print(f"  -> Statut : {status_equateur} (Écart: {diff_equateur:e})")
    print(
        f"  (Note: À l'équateur, thêta s'aligne géométriquement avec -Z. Si {nom_champ}z=0, {nom_champ}_theta doit valoir 0)")

    print("\n" + "=" * 45)

def valider_plan_horizontal_z0(X, Y, F_cartesian, F_theta, F_phi, nom_champ="E"):
    """
    Effectue des tests de cohérence géométrique sur un plan horizontal z=0.

    Parameters:
    -----------
    X, Y : np.ndarray (2D)
        Grilles de coordonnées cartésiennes de forme (ny, nx) ou (nx, ny).
    F_cartesian : np.ndarray (3D)
        Champ cartésien original de forme (ny, nx, 3) ou (nx, ny, 3) -> [Fx, Fy, Fz].
    F_theta : np.ndarray (2D)
        Matrice calculée de la composante theta.
    F_phi : np.ndarray (2D)
        Matrice calculée de la composante phi.
    nom_champ : str
        Nom pour l'affichage ("E" ou "H").
    """
    print(f"=== VALIDATION PLAN HORIZONTAL (Z=0) : CHAMP {nom_champ} ===")

    # Extraction des composantes cartésiennes
    Fx = F_cartesian[..., 0]
    Fy = F_cartesian[..., 1]

    tol = 1e-5

    # -------------------------------------------------------------
    # TEST GLOBALE SUR F_THETA
    # -------------------------------------------------------------
    # À z=0 (theta = pi/2), un champ purement horizontal n'a pas de composante selon theta.
    max_F_theta = np.max(np.abs(F_theta))
    status_theta = "✓ SUCCÈS" if max_F_theta < tol else "❌ ÉCHEC"
    print(f"\n[Test Global] Vérification de {nom_champ}_theta :")
    print(f"  - Valeur maximale absolue de {nom_champ}_theta sur le plan : {max_F_theta:e}")
    print(f"  -> Statut : {status_theta} (Doit être proche de 0 partout car le champ est horizontal à z=0)")

    # -------------------------------------------------------------
    # RECHERCHE DES INDICES DES POINTS CARDINAUX
    # -------------------------------------------------------------
    # Recherche des centres (où X=0 et Y=0)
    idx_center_x = np.argmin(np.abs(X[0, :])) if X.ndim > 1 else np.argmin(np.abs(X))
    idx_center_y = np.argmin(np.abs(Y[:, 0])) if Y.ndim > 1 else np.argmin(np.abs(Y))

    # Extrêmes
    idx_max_x = np.argmax(X[idx_center_y, :])
    idx_min_x = np.argmin(X[idx_center_y, :])
    idx_max_y = np.argmax(Y[:, idx_center_x])
    idx_min_y = np.argmin(Y[:, idx_center_x])

    # -------------------------------------------------------------
    # TEST 1 : À l'Est (X > 0, Y = 0) -> phi = 0
    # -------------------------------------------------------------
    # Le vecteur unitaire phi pointe vers le Nord (+Y). F_phi doit valoir Fy.
    val_Fphi_est = F_phi[idx_center_y, idx_max_x]
    val_Fy_est = Fy[idx_center_y, idx_max_x]
    diff_est = np.abs(val_Fphi_est - val_Fy_est)
    status_est = "✓ SUCCÈS" if diff_est < tol else "❌ ÉCHEC"

    print(f"\n[Test 1] Point EST (X={X[idx_center_y, idx_max_x]:.1f}, Y=0) :")
    print(f"  - {nom_champ}_phi calculé : {val_Fphi_est:e}")
    print(f"  - {nom_champ}_y cartésien : {val_Fy_est:e}")
    print(f"  -> Statut : {status_est} (phi s'aligne avec +Y)")

    # -------------------------------------------------------------
    # TEST 2 : Au Nord (X = 0, Y > 0) -> phi = pi/2
    # -------------------------------------------------------------
    # Le vecteur unitaire phi pointe vers l'Ouest (-X). F_phi doit valoir -Fx.
    val_Fphi_nord = F_phi[idx_max_y, idx_center_x]
    val_Fx_nord = Fx[idx_max_y, idx_center_x]
    diff_nord = np.abs(val_Fphi_nord - (-val_Fx_nord))
    status_nord = "✓ SUCCÈS" if diff_nord < tol else "❌ ÉCHEC"

    print(f"\n[Test 2] Point NORD (X=0, Y={Y[idx_max_y, idx_center_x]:.1f}) :")
    print(f"  - {nom_champ}_phi calculé : {val_Fphi_nord:e}")
    print(f"  - -{nom_champ}_x cartésien : {-val_Fx_nord:e}")
    print(f"  -> Statut : {status_nord} (phi s'aligne avec -X)")

    # -------------------------------------------------------------
    # TEST 3 : À l'Ouest (X < 0, Y = 0) -> phi = pi
    # -------------------------------------------------------------
    # Le vecteur unitaire phi pointe vers le Sud (-Y). F_phi doit valoir -Fy.
    val_Fphi_ouest = F_phi[idx_center_y, idx_min_x]
    val_Fy_ouest = Fy[idx_center_y, idx_min_x]
    diff_ouest = np.abs(val_Fphi_ouest - (-val_Fy_ouest))
    status_ouest = "✓ SUCCÈS" if diff_ouest < tol else "❌ ÉCHEC"

    print(f"\n[Test 3] Point OUEST (X={X[idx_center_y, idx_min_x]:.1f}, Y=0) :")
    print(f"  - {nom_champ}_phi calculé : {val_Fphi_ouest:e}")
    print(f"  - -{nom_champ}_y cartésien : {-val_Fy_ouest:e}")
    print(f"  -> Statut : {status_ouest} (phi s'aligne avec -Y)")

    # -------------------------------------------------------------
    # TEST 4 : Au Sud (X = 0, Y < 0) -> phi = 3*pi/2
    # -------------------------------------------------------------
    # Le vecteur unitaire phi pointe vers l'Est (+X). F_phi doit valoir Fx.
    val_Fphi_sud = F_phi[idx_min_y, idx_center_x]
    val_Fx_sud = Fx[idx_min_y, idx_center_x]
    diff_sud = np.abs(val_Fphi_sud - val_Fx_sud)
    status_sud = "✓ SUCCÈS" if diff_sud < tol else "❌ ÉCHEC"

    print(f"\n[Test 4] Point SUD (X=0, Y={Y[idx_min_y, idx_center_x]:.1f}) :")
    print(f"  - {nom_champ}_phi calculé : {val_Fphi_sud:e}")
    print(f"  - {nom_champ}_x cartésien : {val_Fx_sud:e}")
    print(f"  -> Statut : {status_sud} (phi s'aligne avec +X)")

    print("\n" + "=" * 50)
