from typing import Tuple
import pandas as pd
import matplotlib.pyplot as plt

from delta_inversion import *


def test_delta_inv(delta_fwd: float):
    vol_rng = np.linspace(0.001, 1.0, 50)
    for i, vol in enumerate(vol_rng):
        z_call = z_prem_adj(vol, delta_fwd, OptionType.CALL)
        z_put = z_prem_adj(vol, -delta_fwd, OptionType.PUT)
        assert abs(delta_prem_adj(vol, z_call, OptionType.CALL) - delta_fwd) < 1e-8
        assert abs(delta_prem_adj(vol, z_put, OptionType.PUT) + delta_fwd) < 1e-8


def test_r_inv():
    ys = [0.1, 0.5, np.sqrt(np.pi / 2), 2.0, 10.0]
    for y in ys:
        x = R_inv(y)
        assert np.isfinite(x) and abs(R(x) - y) < 1e-8


def test_delta_prem_adj():
    total_vols = np.geomspace(1e-3, 2.0, 50)
    deltas = np.geomspace(1e-4, 0.5, 30)
    for total_vol in total_vols:
        for d in deltas:
            z_call = z_prem_adj(total_vol, d, OptionType.CALL)
            if np.isfinite(z_call):
                assert abs(delta_prem_adj(total_vol, z_call, OptionType.CALL) - d) < 1e-8
            z_put = z_prem_adj(total_vol, -d, OptionType.PUT)
            if np.isfinite(z_put):
                assert abs(delta_prem_adj(total_vol, z_put, OptionType.PUT) + d) < 1e-8


def test_delta_prem_adj_with_tol(tol: float) -> Tuple[float, float]:
    ttxs = np.atleast_1d([1.0 / 365, 1 / 12, 1.0, 10.0])
    vols = np.linspace(0.01, 5.0, 10)
    deltas = np.geomspace(1e-4, 0.5, 30)

    mse_calls = 0.0
    mse_puts = 0.0
    nb_calls = 0
    nb_puts = 0
    for vol in vols:
        for ttx in ttxs:
            total_vol = vol * np.sqrt(ttx)
            for d in deltas:
                z_call = z_prem_adj(total_vol, d, OptionType.CALL, tol=tol)
                z_call_bm = z_prem_adj(total_vol, d, OptionType.CALL, tol=1e-12)
                if np.isfinite(z_call):
                    assert abs(delta_prem_adj(total_vol, z_call, OptionType.CALL) - d) < tol
                    mse_calls += (z_call - z_call_bm) ** 2
                    nb_calls += 1
                z_put = z_prem_adj(total_vol, -d, OptionType.PUT, tol=tol)
                z_put_bm = z_prem_adj(total_vol, -d, OptionType.PUT, tol=1e-12)
                if np.isfinite(z_put):
                    assert abs(delta_prem_adj(total_vol, z_put, OptionType.PUT) + d) < tol
                    mse_puts += (z_put - z_put_bm) ** 2
                    nb_puts += 1

    rmse_calls = np.sqrt(mse_calls / nb_calls)
    rmse_puts = np.sqrt(mse_puts / nb_puts)

    return (rmse_puts, rmse_calls)


def test_delta_prem_rmse():
    rows = []
    for tol in [1e-6, 1e-8, 1e-9]:
        rmse_puts, rmse_calls = test_delta_prem_adj_with_tol(tol=tol)
        rows.append({'Tolerance': tol, 'RMSE puts': rmse_puts, 'RMSE calls': rmse_calls})

    return rows


def test_delta_max(tol: float):
    ttxs = np.atleast_1d([1.0 / 365, 1 / 12, 1.0, 10.0])
    vols = np.linspace(0.01, 5.0, 10)
    alphas = [0.99, 0.999, 0.9999]
    mse = 0
    for vol in vols:
        for ttx in ttxs:
            total_vol = vol * np.sqrt(ttx)
            fo = R_inv(1.0 / total_vol)
            z_max = total_vol * (fo - 0.5 * total_vol)
            delta_max = delta_prem_adj(total_vol, z_max, OptionType.CALL)
            for alpha in alphas:
                delta = alpha * delta_max
                z = z_prem_adj(total_vol, delta, OptionType.CALL, tol=tol)
                z_bm = z_prem_adj(total_vol, delta, OptionType.CALL, tol=1e-12)
                assert np.isfinite(z)
                assert abs(delta_prem_adj(total_vol, z, OptionType.CALL) - delta) < tol
                mse += (z - z_bm) ** 2
    nb_calls = vols.size * ttxs.size * len(alphas)
    rmse = np.sqrt(mse / nb_calls)
    return rmse


def test_delta_max_neighbourhood():
    rows = []
    for tol in [1e-6, 1e-8, 1e-9]:
        rmse = test_delta_max(tol=tol)
        rows.append({'Tolerance': tol, 'RMSE': rmse})

    return rows


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

    fig, axs = plt.subplots(1, 2, figsize=(10, 4.0), tight_layout=True)

    axs[0].plot(z_rng, delta_calls)
    axs[0].axhline(y=delta_max + 2e-3, color='red')  # shift slightly up for clarity
    axs[0].axvline(x=z_max, color='green')
    axs[0].text(z_rng[0], delta_max, r'$\Delta_{max}$', color='red',
                va='bottom', ha='left', fontsize=10)
    axs[0].text(z_max, axs[0].get_ylim()[0] if False else min(delta_calls),
                r'$z_{max}$', color='green', va='bottom', ha='right', fontsize=10)

    axs[0].set_xlabel(f"$z$")
    axs[0].set_ylabel(f"$\Delta^{{pa,F}}$")
    axs[0].set_title(f"Premium-adjusted forward delta $\\Delta^{{pa,F}}$\n as a function of log-moneyness $z$",
                     fontsize=10, color="darkblue")

    axs[1].plot(delta_rng, z_pa_rng)
    axs[1].set_xlabel(f"$\Delta^{{pa,F}}$")
    axs[1].set_ylabel(f"$z$")
    axs[1].set_title(f"Log-moneyness $z$ as a function of \n premium-adjusted forward delta $\\Delta^{{pa,F}}$",
                     fontsize=10, color="darkblue")

    plt.show()
    fig.savefig("figures/prem_adj_delta.png", dpi=150, bbox_inches='tight')



def plot_W_deriv(delta_fwd: float):
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.0), tight_layout=True)
    vols = np.linspace(0.25, 0.5, 101)
    T = 6
    total_vols = vols * np.sqrt(T)
    W_r_derivs = np.array([W_deriv(total_vol, delta_fwd, True) for total_vol in total_vols])
    W_l_derivs = np.array([W_deriv(total_vol, delta_fwd, False) for total_vol in total_vols])
    axs[0].plot(vols, W_r_derivs)
    axs[1].plot(vols, W_l_derivs)

    axs[0].set_xlabel(f"$\\sigma_{{\mathrm{{str}}}}$")
    # axs[0].set_ylabel(f"$\\Delta^{{pa,F}}$")
    axs[0].set_title(f"Derivative $W'(\sigma_{{\mathrm{{str}}}})$ as a function of $\sigma_{{\mathrm{{str}}}}$,\n"
                     f"assuming smaller call root $z_a$",
                     fontsize=10, color="darkblue")

    axs[0].set_xlabel(f"$\sigma_{{\mathrm{{str}}}}$")
    # axs[1].set_ylabel(f"$z$")
    axs[1].set_title(f"Derivative $W'(\sigma_{{\mathrm{{str}}}})$ as a function of $\sigma_{{\mathrm{{str}}}}$,\n"
                     f"assuming larger call root $z_b$",
                     fontsize=10, color="darkblue")

    plt.show()
    fig.savefig("figures/strangle_vega.png", dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    print(f"Running 'test_r_inv':")
    test_r_inv()

    print(f"Running 'test_delta_inv':")
    test_delta_inv(delta_fwd=0.01)

    print(f"Running 'test_delta_prem_adj':")
    test_delta_prem_adj()

    print(f"Running 'plot_prem_adj_delta':")
    plot_prem_adj_delta()

    print(f"Running 'plot_W_deriv':")
    delta_fwd = 0.25
    plot_W_deriv(delta_fwd)

    print(f"Running 'test_delta_prem_rmse':")
    df1 = pd.DataFrame(test_delta_prem_rmse()).to_string(formatters={
        'Tolerance': '{:.2e}'.format,
        'RMSE puts': '{:.2e}'.format,
        'RMSE calls': '{:.2e}'.format
    })
    print(df1)

    print(f"Running 'test_delta_max_neighbourhood':")
    df2 = pd.DataFrame(test_delta_max_neighbourhood()).to_string(formatters={
        'Tolerance': '{:.2e}'.format,
        'RMSE': '{:.2e}'.format
    })
    print(df2)
