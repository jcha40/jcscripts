/**
    * @param {Array} controlpoints - Array of control points
    * @return {Function} - Function that takes a parameter t and returns a point on the bezier curve [x, y]
 */
let bezier = function(...controlpoints) {
    let n = controlpoints.length - 1,
        coefs = Array(n + 1).fill(1);

    for (let i = 0; i < n; i++) {
        coefs[i + 1] = coefs[i] * (n - i) / (i + 1)
    };

    return function(t) {
        let x = 0,
            y = 0;
        for (let i = 0; i <= n; i++) {
            let b = coefs[i] * Math.pow(1 - t, n - i) * Math.pow(t, i);
            x += b * controlpoints[i][0];
            y += b * controlpoints[i][1];
        }
        return [x, y]
    }
};

/**
    * @param {Array} controlpoints - Array of control points
    * @return {Function} - Function that takes a parameter t and returns the derivative of the bezier curve at that point [dx, dy]
 */
let bezier_derivative = function(...controlpoints) {
    let bezier1 = bezier(...controlpoints.slice(1)),
        bezier2 = bezier(...controlpoints.slice(0, -1));

    return function(t) {
        let b1 = bezier1(t),
            b2 = bezier2(t);
        return [b1[0] - b2[0], b1[1] - b2[1]]
    }
}