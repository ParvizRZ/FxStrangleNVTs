import numpy as np
from delta_inversion import *


def test_delta_inv(delta_fwd: float):
    vol_rng = np.linspace(0.001, 1.0, 50)
    for i, vol in enumerate(vol_rng):
        z_call = z_prem_adj(vol, delta_fwd, OptionType.CALL)
        z_put = z_prem_adj(vol, -delta_fwd, OptionType.PUT)
        assert abs(delta_prem_adj(vol, z_call, OptionType.CALL) - delta_fwd) < tol
        assert abs(delta_prem_adj(vol, z_put, OptionType.PUT) + delta_fwd) < tol


def test_r_inv():
    ys = [0.1, 0.5, np.sqrt(np.pi / 2), 2.0, 10.0]
    for y in ys:
        x = R_inv(y)
        assert np.isfinite(x) and abs(R(x) - y) < tol


def test_delta_prem_adj():
    total_vols = np.geomspace(1e-3, 2.0, 50)
    deltas = np.geomspace(1e-4, 0.5, 30)
    for total_vol in total_vols:
        for d in deltas:
            z_call = z_prem_adj(total_vol, d, OptionType.CALL)
            if np.isfinite(z_call):
                assert abs(delta_prem_adj(total_vol, z_call, OptionType.CALL) - d) < tol
            z_put = z_prem_adj(total_vol, -d, OptionType.PUT)
            if np.isfinite(z_put):
                assert abs(delta_prem_adj(total_vol, z_put, OptionType.PUT) + d) < tol


def plot_prem_adj_delta():
    vol = 0.15
    ttx = 2
    total_vol = vol * np.sqrt(ttx)
    z_rng = np.linspace(-2.0, 1.5, 501) * total_vol
    delta_calls = [delta_prem_adj(total_vol, z, OptionType.CALL) for z in z_rng]

    fo = R_inv(1.0 / total_vol)
    z_max = total_vol * (fo - 0.5 * total_vol)
    delta_max = delta_prem_adj(total_vol, z_max, OptionType.CALL)

    delta_rng = np.linspace(1e-4, delta_max - 1e-4, 501)
    z_pa_rng = [z_prem_adj(total_vol, delta, OptionType.CALL) for delta in delta_rng]

    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.0), tight_layout=True)

    axs[0].plot(z_rng, delta_calls)
    axs[0].axhline(y=delta_max+2e-3, color='red')  # shift slightly up for clarity
    axs[0].axvline(x=z_max, color='green')
    axs[0].text(z_rng[0], delta_max, r'$\Delta_{max}$', color='red',
                va='bottom', ha='left', fontsize=10)
    axs[0].text(z_max, axs[0].get_ylim()[0] if False else min(delta_calls),
                r'$z_{max}$', color='green', va='bottom', ha='right', fontsize=10)


    axs[0].set_xlabel(f"$z$")
    axs[0].set_ylabel(f"$\Delta^{{pa,F}}$")
    axs[0].set_title(f"Premium-adjusted forward delta $\\Delta^{{pa,F}}$\n as a function of log-moneyness $z$", fontsize=10, color="darkblue")

    axs[1].plot(delta_rng, z_pa_rng)
    axs[1].set_xlabel(f"$\Delta^{{pa,F}}$")
    axs[1].set_ylabel(f"$z$")
    axs[1].set_title(f"Log-moneyness $z$ as a function of \n premium-adjusted forward delta $\\Delta^{{pa,F}}$", fontsize=10, color="darkblue")

    # plt.savefig("figures/prem_adj_delta.png", dpi=150, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    print(f"Running 'test_r_inv':")
    test_r_inv()

    print(f"Running 'test_delta_inv':")
    test_delta_inv(delta_fwd=0.01)

    print(f"Running 'test_delta_prem_adj':")
    test_delta_prem_adj()

    print(f"Running 'plot_prem_adj_delta':")
    plot_prem_adj_delta()
