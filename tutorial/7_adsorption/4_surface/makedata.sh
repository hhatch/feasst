awk '{if ($1 == "O") {q=-0.75; type=1} else {q=1.5; type=2}; print NR, "1", type, q, $2, $3, $4, "0 0 0"}' ttt > tt
