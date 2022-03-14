/**
 * Piecewise cubic Hermite interpolating polynomial
 * @param   {Array} u   Array of values to interpolate on; must be sorted
 * @param   {Array} x   Array containing input domain values; must be sorted
 * @param   {Array} y   Array containing input range values
 * @return  {Array}     Array containing interpolated range values at domain points in u
 */
function pchip(u, x, y) {
    let h = [];
    for (let i in x) {
        h.push(x[i] - x[i - 1]);
    }

    let n = 0,
        y_hat = [];
    for (let i in u) {
        if (u[i] < x[0] || u[i] > x[x.length - 1]) {
            y_hat.push(NaN);
            continue
        }
        while (x[n] <= u[i]) {
            n++;
        }
        let p0 = y[n - 1],
            p1 = y[n],
            d0 = (p0 - y[n - 2]) / h[n - 1],
            d1 = (p1 - p0) / h[n],
            d2 = (y[n + 1] - p1) / h[n + 1],
            fp0 = n === 1 ? d1 : (h[n - 1] * d0 + h[n] * d1) / (h[n - 1] + h[n]),
            fp1 = n === (x.length - 1) ? d1 : (h[n] * d1 + h[n + 1] * d2) / (h[n] + h[n + 1]),
            v = (u[i] - x[n - 1]) / (x[n] - x[n - 1]);
        y_hat.push((p0 * ((2 * v * v * v) - (3 * v * v) + 1)) +
            (p1 * ((-2 * v * v * v) + (3 * v * v))) +
            (fp0 * ((v * v * v) - (2 * v * v) + v)) +
            (fp1 * ((v * v * v) - (v * v))));
    }
    return y_hat
}
