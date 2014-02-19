(ns linear-systems.core-test
  (:require [clojure.test :refer :all]
            [linear-systems.core :as ls :refer :all]))

(deftest basic-inverse-matrix-test
  ; This test includes functions
	;  - lu-factorization-sparse
	;  - pivot
	;  - lu-forward-substitution-sparse
	;  - backward-substitution
	;  - inverse-matrix
	;  - diag
	;  - matrix-prod
  (testing "Basic inverse-matrix."
    (let [A  [[1 2 3 4][6 5 4 1][9 8 7 1][9 2 1 3]]
          invA (ls/inverse-matrix A)]
           (is (= (diag (repeat (count A) 1)) (matrix-prod A invA))))))
