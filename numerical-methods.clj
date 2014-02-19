;(ns
;  ^{:author "Antti Rämö",
;    :doc ""}
;  numerical-methods)

;(def A [[1 2 3 4][6 5 4 1][9 8 7 1][9 2 1 3]])
;(def A [[1 2 3][4 5 1][3 2 1]]) (def b [8 5 2])
;(def A [[1.0E-32 1][1 1]]) (def b [1 2])
;(def A [[0.6 1.52 3.5][2 4 1][1 2.8 1]])
;(def A [[1 2 3][0 4 5][1 0 6]])
;(def A [[4 4 5 5][1 2 3 4][3 2 1 1][8 2 1 7]])
;(def A [[8 2 1 7][4 4 5 5][3 2 1 1][1 2 3 4]])
;(def A2 [[4 6 2][6 18 -1.5][2 -1.5 4.25]])
;(def A [[1 2 3 8][8 7 5 2][9 6 2 1][7 5 3 1]])
(defn pivot1 [A k]
  (reduce 
    #(if (< (second %2) (nth (first %1) k))
      %1 
      (cons 
        (nth %1 (first %2))
        (concat (take (first %2) %1) (drop (inc (first %2)) %1))))
    A (map vector (range) (map #(nth % k) A))))
;(println (pivot1 A 1))

(defn diag [A] 
  ; (diag A) Returns a diagonal element of a matrix
  ; (diag v) Returns a diagonal matrix with diagonal units v.
  (if (coll? (first A))
    (map nth A (range))
    (map #(concat (repeat %2 0) [%1] (repeat (- (count A) (inc %2)) 0)) A (range))))

(defn matrix-vector-prod [m v] 
  (map #(apply + (map * % v)) m))

(defn matrix-prod 
  "Calculates matrix production A*B. Function takes two matrices
  A and Bt, where Bt is transpose of B. When A is size n*k and
  Bt is size m*k returned matrix is size n*m. Matrices are 
  represented by vector of row vectors."
  [A Bt]
  (mapv (fn [ai] (mapv (fn [bj] (apply + (map * ai bj))) Bt)) A))

(defn triangular-to-full
  "Converts sparse triangular matrix A to full matrix where empty
  elements are filled with zeros. Matrices are represented by 
  vector of row vectors."
  [A]
  (let [n (count A)]
    (if (< (count (first A)) n)
      (map #(into (vec %) (repeat (- n (count %)) 0)) A)    ; Lower triangular matrix.
      (map #(into (vec (repeat (- n (count %) 0)) %) A))))) ; Upper triangular matrix.

(defn pivot
  ([A] (pivot A (range (count A)) 0))
  ([A k] (pivot A (range (count A)) k))
  ([A p k]
    ; Find the row with the greatest value in the column k.
    (let [i (->> 
              ;(map first A)  ; Take first column.
              (map #(nth % k) A)  ; Take column k.
              (map #(* % %))      ; Square piecewise
              (zipmap (range))    ; Add row index.
              (apply max-key val) ; Take the key-value-pair with greatest value.
              (key))              ; Extract the key.
          r (into [i] (remove (partial = i) (range (count A))))] ; New row order.
;          r (assoc (vec (range (count A))) 0 i i 0)] ; New row order.
    ; Rearrange the data
      {:A (replace (vec A) r)
       :p (replace (vec p) r)
       :r r})))
;(pivot A 1)

;(defn backward-substitution [U c] 
;; Solves a linear system 'Ux = c' by calculating the backward substitution.
;    ; Take upper triangular matrix without zeros and cons vector c.
;    (let [U (map
;              #(cons %1 (take-last %3 %2)) 
;              c U (iterate dec (count U)))]
;      (reduce
;        #(cons 
;          (/ (- (first %2) (apply + (map * %1 (nnext %2)))) (second %2))
;          %1)
;        '() (reverse U))))
;(println (backward-substitution [[2 4 1][0 0.8 0.5][0 0 3]] [-3.8 -0.4 9.6]))

(defn backward-substitution-sparse [U c] 
; Solves a linear system 'Ux = c' by calculating the backward substitution.
  (if (map? U)
    ; If U is returned from LU-factorization (returns map: {:L L :U U :p p}).
    (backward-substitution-sparse (U :U) c)
    (if (coll? (first c))
      ; If c is a matrix calculate each column separately.
      (map (partial backward-substitution-sparse U) c)
      ; Take upper triangular matrix without zeros and cons vector c.
      (let [U (map cons c U)]
        (reduce
          #(cons 
            (/ (- (first %2) (apply + (map * %1 (nnext %2)))) (second %2))
            %1)
          '() (reverse U))))))


(defn gaussian-elimination
  ([A b] 
    (gaussian-elimination (map #(conj (vec %1) %2) A b)))
  ([A]
    (letfn [(elimination-step [A]
              (map 
                (fn [ai-jj]
                  (map 
                    #(- %1 (/ (* %2 (first ai-jj)) (first (first A))))
                    (rest ai-jj) (rest (first A)))) 
                (rest A)))
            (sub-ge [A k]
              (when-first [a1 A]
                (cons (concat (repeat k 0) a1)
                      (sub-ge (elimination-step A) (inc k)))))]
    (sub-ge A 0))))
;(println (gaussian-elimination A b))
;(println "gauss A2: " (gaussian-elimination A2))

(defn gaussian-elimination-with-pivoting 
  ([A b]
    (gaussian-elimination-with-pivoting (map #(conj (vec %1) %2) A b)))
  ([A]
    (letfn [(elimination-step [A]
              (map 
                (fn [ai]
                  (map 
                    #(- %1 (/ (* %2 (first ai)) (first (first A))))
                    (rest ai) (rest (first A)))) 
                (rest A)))
            (sub-ge [A k]
              (when-first [a1 A]
                (cons (concat (repeat k 0) a1)
                      (sub-ge (elimination-step ((pivot A) :A)) (inc k)))))]
    (sub-ge A 0))))
;(println (gaussian-elimination-with-pivoting A b))

;(defn lu-factorization [A]
;  (letfn [(elimination-step [A]
;            (map 
;              (fn [ai]
;                (let [m (/ (first ai) (first (first A)))]
;                  [(map #(- %1 (* %2 m)) (rest ai) (rest (first A))) m]))
;              (rest A)))
;          (sub-ge [A k p]
;            (when-first [a1 A]
;              (let [Ap (pivot A p 0)
;                    Am (elimination-step (Ap :A))]
;                (cons [(concat (repeat k 0) [1] (reverse (map second Am))) ; Li (transpose)
;                       (concat (repeat k 0) (first (Ap :A)))     ; return Ui
;                       (first (Ap :p))]                                ; pi
;                      (sub-ge (map first Am) (inc k) (rest (Ap :p)))))))]
;  (apply map vector (sub-ge A 0 (range (count A))))))
;(def LUp (lu-factorization A))
;(println "LUp = " LUp)



(defn lu-factorization-sparse [A]
  (letfn [(elimination-step [[A1 & Ak]]
            (reduce 
              (fn [AL [ai1 & aik]]
                (let [l (/ ai1 (first A1))
                      ai (mapv #(- %1 (* %2 l)) aik (rest A1))]
                  {:A (conj (AL :A) ai) :L (conj (AL :L) l)}))
              {:A [] :L []} Ak))
          (rearrange-L [L p]
            (let [pn (rest (reductions
                      (fn [[p1 & pk] i] (mapv #(if (< p1 %) (dec %) %) pk))
                      p p))]
              (mapv replace L pn)))]
    (loop [A A i (range (count A)) L [] U [] p []]
      (if (first A)
        (let [Ap (pivot A i 0)
              Al (elimination-step (Ap :A))
              lj (Al :L)
              ui (first (Ap :A))
              pi (first (Ap :p))]
          (recur (Al :A) (rest (Ap :p)) (conj L lj) (conj U ui) (conj p pi)))
      (zipmap [:L :U :p] [(rearrange-L L p) U p])))))
;(def LUp (lu-factorization-sparseb A))
;(def LUp2 (lu-factorization-sparse2 A2))
;(println "LUp2 = " LUp)
;(println "LUp2b = " LUp2)


(defn lu-factorization-sparse2 [A]
  (letfn [(elimination-step [[[a11 & a1k] & Ak]]
            (->>
              (map (fn [[ak1 & aki]]
                      (let [l (/ ak1 a11)
                            ak (mapv #(- %1 (* %2 l)) aki a1k)]
                        [ak l]))       ; Return values for each row
                    Ak)
              (apply map vector)       ; Separate A and L parts
              (if (not Ak) [[][]])))   ; If Ak is empty return empty structure
          (rearrange-L [L p]
            (let [pn (rest (reductions
                      (fn [[p1 & pk] i] (mapv #(if (< p1 %) (dec %) %) pk))
                      p p))]
              (mapv replace L pn)))]
    (loop [A A i (range (count A)) L [] U [] p []]
      (if (first A)
        (let [Ap (pivot A i 0)
              Al (elimination-step (Ap :A))
              lj (second Al)
              ui (first (Ap :A))
              pi (first (Ap :p))]
          (recur (first Al) (rest (Ap :p)) (conj L lj) (conj U ui) (conj p pi)))
      (zipmap [:L :U :p] [(rearrange-L L p) U p])))))
;(def LUp (lu-factorization-sparseb A))
;(def LUp2 (lu-factorization-sparse2 A2))
;(println "LUp2 = " LUp)
;(println "LUp2b = " LUp2)


;(defn lu-factorization-sparseb [A]
;  (letfn [; Recursive step based on gaussian elimination.
;          (elimination-step [[[a11 & a1k] & Ak]]
;            (reduce 
;              (fn [AL [ai1 & aik]]
;                (let [l (/ ai1 a11)
;                      ai (mapv #(- %1 (* %2 l)) aik a1k)]
;                  {:A (conj (AL :A) ai) :L (conj (AL :L) l)}))
;              {:A [] :L []} Ak))
;          ; Rearrange each column of L matrix by pivoting index.
;          (rearrange-L [L p]
;              (mapv 
;                #(replace %1 (into (vec (range %2)) (map (partial + %2) p)))
;                L (range (dec (count L)) -1 -1)))]
;    (loop [A A i (range (count A)) L [] U [] p []]
;      (if (first A)
;        (let [Ap (pivot A i 0)              ; Take the pivot
;              Al (elimination-step (Ap :A)) ; Calculate the elimination
;              L (rearrange-L L (Ap :r))     ; Sort the values on L
;              lj (Al :L)                    ; Take new col on L
;              ui (first (Ap :A))            ; Take new row on U
;              pi (first (Ap :p))]           ; Take new value on p
;          (recur (Al :A) (rest (Ap :p)) (conj L lj) (conj U ui) (conj p pi)))
;      (zipmap [:L :U :p] [L U p])))))
;(def LUp (lu-factorization-sparse A))
;;(def LUp2 (lu-factorization-sparse2 A2))
;(println "LUp = " LUp)
;;(println "LUp2b = " LUp2)


;(defn lu-factorization-sparseb2 [A]
;  (letfn [; Recursive step based on gaussian elimination.
;          (elimination-step [[[a11 & a1k] & Ak]]
;            (->>
;              (map (fn [[ak1 & aki]]
;                      (let [l (/ ak1 a11)
;                            ak (mapv #(- %1 (* %2 l)) aki a1k)]
;                        [ak l]))       ; Return values for each row
;                    Ak)
;              (apply map vector)       ; Separate A and L parts
;              (if (not Ak) [[][]])))   ; If Ak is empty return empty structure
;          (rearrange-L [L p]
;              (mapv 
;                #(replace %1 (into (vec (range %2)) (map (partial + %2) p)))
;                L (range (dec (count L)) -1 -1)))]
;    (loop [A A i (range (count A)) L [] U [] p []]
;      (if (first A)
;        (let [Ap (pivot A i 0)              ; Take the pivot
;              Al (elimination-step (Ap :A)) ; Calculate the elimination
;              L (rearrange-L L (Ap :r))     ; Sort the values on L
;              lj (second Al)                ; Take new col on L
;              ui (first (Ap :A))            ; Take new row on U
;              pi (first (Ap :p))]           ; Take new value on p
;          (recur (first Al) (rest (Ap :p)) (conj L lj) (conj U ui) (conj p pi)))
;      (zipmap [:L :U :p] [L U p])))))
;(def LUp (lu-factorization-sparse2 A))
;;(def LUp2 (lu-factorization-sparse2 A2))
;(println "LUp = " LUp)
;;(println "LUp2b = " LUp2)






(defn lu-forward-substitution-sparse
  ([LUP b] (lu-forward-substitution-sparse (LUP :L) (LUP :p) b))
  ([L p b]
    (if (coll? (first b))
      ; If b is a matrix (transposed) calculate each column separately.
      (map (partial lu-forward-substitution-sparse L p) b)
      (let [b (map (partial nth b) p)]
        (reduce
          (fn [b k]
            (concat
              (take k b)
              (map 
                #(- %1 (* %2 (nth b (dec k)))) 
                (drop k b) (nth L (dec k)))))
          b (range 1 (count b)))))))

;(def y (lu-forward-substitution-sparse LUp [[8.3 -3.8 -2.3][5 6 2][1 8 1]]))
;(println "y: " y)

;(println "x: " (backward-substitution-sparse LUp y))

(defn inverse-matrix [A]
  (let [LUP (lu-factorization-sparse A)]   ; Calculate LU factorization
    (->>
      (diag (repeat (count A) 1))          ; Create unit matrix
      (lu-forward-substitution-sparse LUP) ; Evaluate forward substition
      (backward-substitution-sparse LUP)   ; Evaluate backward substition
      (apply map vector))))                ; Take transpose
;(println (inverse-matrix A))


(defn cholesky-decomposition
  "Calculates Cholesky decomposition by The Cholesky–Banachiewicz
  algorithm. Function returns lower triagular matrix C, so that 
  A = CC*, where input A is symmetric positive-definite matrix. 
  Matrices are represented by vector of row vectors."
  [A]
  (reduce
    ; Calculate C row-by-row
    (fn [C ai]
      (let [lower-units ; Calculate lower triangular matrix elements of row i.
              (reduce 
                (fn [ci cj] ; Calculate row ci with row cj and row ai.
                  (conj ci (/ (- (nth ai (count ci))
                                 (apply + (map * ci cj)))
                              (last cj))))
                [] C)
            diag-unit ; Calculate diagonal element.
              (Math/sqrt (- (last ai) ; Calculate cii with aii and row ci.
                            (apply + (map #(* % %) lower-units))))]
        (conj C (conj lower-units diag-unit))))
    [] (map take (range 1 (inc (count A))) A)))
;(cholesky-decomposition [[25 15 -5][15 18 0][-5 0 11]])
(def C 
  (triangular-to-full
    (cholesky-decomposition
      [[18 22 54 42][22 70 86 62][54 86 174 134][42 62 134 106]])))
(matrix-prod C C)
