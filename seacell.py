import os
import pickle

import scanpy as sc
import SEACells

data_path = "/home/jiehoonk/mnt/dataset/GBmap/"
out_path = "extended/"


adata = sc.read_h5ad(os.path.join(out_path, "total.h5ad"))
adata.X = adata.layers["scaled"]

n_SEACells = int(adata.shape[0] / 100)
build_kernel_on = "X_umap"
n_waypoint_eigs = 10

model = SEACells.core.SEACells(
    adata,
    build_kernel_on=build_kernel_on,
    n_SEACells=n_SEACells,
    n_waypoint_eigs=n_waypoint_eigs,
    convergence_epsilon=1e-5,
    use_sparse=True,
)

model.construct_kernel_matrix()
M = model.kernel_matrix

print("model construction finished")
with open(os.path.join(out_path, "model/model_150.pkl"), "wb") as f:
    pickle.dump(model, f)


with open(os.path.join(out_path, "model/M.pkl"), "wb") as f:
    pickle.dump(M, f)


model.initialize_archetypes()
print("model initialization finished")

with open(os.path.join(out_path, "model/model_init_150.pkl"), "wb") as f:
    pickle.dump(model, f)


adata.write(os.path.join(out_path, "seacell.h5ad"))

model.fit(min_iter=10, max_iter=50) # 여기서 안됩니다 ㅜ
print(f"Ran for {len(model.RSS_iters)} iterations")
print("model running fished")

with open(os.path.join(out_path, "model/model_train_150.pkl"), "wb") as f:
    pickle.dump(model, f)

adata.write(os.path.join(out_path, "seacell.h5ad"))
