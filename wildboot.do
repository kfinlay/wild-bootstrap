// certification script for wildboot.ado
cscript wildboot adofile wildboot

use xmpl
mycmd x1 x2
assert abs(r(z)-2.5)<=1e-15

keep if x3==2
mycmd x1 x2
local hold = r(z)

use xmpl, clear
mycmd x1 x2 if x3==2
assert r(z) == `hold'

rcof "noisily mycmd x1" == 102       /* too few variables specified */
rcof "noisily mycmd x1 x2 x3" == 103 /* too many variables specified */

// Writing a good certification script
// 
// 
// The purpose of a certification script is to
// 
// 
//     1.  test that the command produces the right answers in a few cases
//         where the solution is known;
// 
// 
//     2.  establish that the command works as it should under extreme
//         conditions, such as R^2=1 regressions;
// 
// 
//     3.  verify that the command responds to mistakes users are likely to
//         make in a dignified manner.
// 
// 
// Certification scripts are written for two reasons:
// 
// 
//     1.  To establish on day one (the day the command is written) that the
//         command works.
// 
// 
//     2.  To allow future changes to be made to the command with some
//         assurance that the changes really are improvements.  (One simply
//         reruns the certification script.)
