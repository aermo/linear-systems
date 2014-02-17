;(def A [[1 2 3 4][6 5 4 1][9 8 7 1][9 2 1 3]])
;(def A [[1 2 3][4 5 1][3 2 1]]) (def b [8 5 2])
;(def A [[1.0E-32 1][1 1]]) (def b [1 2])
(def A [[0.6 1.52 3.5][2 4 1][1 2.8 1]])
(def A2 [[4 6 2][6 18 -1.5][2 -1.5 4.25]])
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
    (map (partial nth) A (range))
    (map #(concat (repeat %2 0) [%1] (repeat (- (count A) (inc %2)) 0)) A (range))))

(defn matrix-vector-prod [m v] 
  (map #(apply + (map * % v)) m))

(defn matrix-prod [A B]
  (map (fn [ai] (map (fn [bj] (apply + (map * ai bj))) (apply map list B))) A))

(defn pivot [A k]
  ; Find the row with the greatest value in the column k.
  (let [p (->> 
            (map first A)  ; Take column k.
            ;(map #(nth % k) A)  ; Take column k.
            (zipmap (range))    ; Add row index.
            (apply max-key val) ; Take the key-value-pair with greatest value.
            (key))]             ; Extract the key.
  ; Move the row on the top
    {:A (concat [(nth A p)] (take p A) (drop (inc p) A)) :p p}))
;(println (pivot A 1))

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
  (if (coll? (first c))
    ; If c is a matrix (transposed) calculate each column separately.
    (map (partial backward-substitution-sparse U) c) ; Then
    ; Take upper triangular matrix without zeros and cons vector c.
    (let [U (map cons c U)] ; Else
      (reduce
        #(cons 
          (/ (- (first %2) (apply + (map * %1 (nnext %2)))) (second %2))
          %1)
        '() (reverse U)))))


(defn backward-substitution1 [U c] 
; Solves a linear system 'Ux = c' by calculating the backward substitution.
    ; Take upper triangular matrix without zeros and cons vector c.
    (let [U (map
              #(cons %1 (take %3 (reverse %2))) 
              c U (iterate dec (count U)))]
      (reduce
        #(cons 
          (/ (- (first %2) (apply + (map * %1 (rest %2)))) (last %2))
          %1)
        '() (reverse U))))
;(println (backward-substitution [[2 4 1][0 0.8 0.5][0 0 3]] [-3.8 -0.4 9.6]))

(defn backward-substitution2 [U c] 
; Backward substitution with a loop structure.
  (loop [k (dec (count U)) y '()]
    (if (< k 0)
      y
      (recur
        (dec k)
        (cons 
          (/ 
            (- 
              (nth c k) 
              (apply + (map * (take-last (- (dec (count U)) k) (nth U k)) y)))
            (nth (nth U k) k)) y)))))
;(println (backward-substitution2 [[2 4 1][0 0.8 0.5][0 0 3]] [-3.8 -0.4 9.6]))


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
(println "gauss A2: " (gaussian-elimination A2))

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
                      (sub-ge (elimination-step ((pivot A 0) :A)) (inc k)))))]
    (sub-ge A 0))))
;(println (gaussian-elimination-with-pivoting A b))

(defn lu-factorization [A]
  (letfn [(elimination-step [A]
            (map 
              (fn [ai]
                (let [m (/ (first ai) (first (first A)))]
                  [(map #(- %1 (* %2 m)) (rest ai) (rest (first A))) m]))
              (rest A)))
          (sub-ge [A k p]
            (when-first [a1 A]
              (let [Ap (pivot A 0)
                    Am (elimination-step (Ap :A))]
                (cons [(concat (repeat k 0) [1] (reverse (map second Am))) ; Li (transpose)
                       (concat (repeat k 0) (first (Ap :A)))     ; return Ui
                       (nth p (Ap :p))]                                ; pi
                      (sub-ge (map first Am) (inc k) (remove (partial = (Ap :p)) p))))))]
  (apply map vector (sub-ge A 0 (range (count A))))))
(def LUp (lu-factorization A))
(println "LUp = " LUp)

(defn lu-factorization-sparse [A]
  (letfn [(elimination-step [A]
            (reduce 
              (fn [AL ai]
                (let
                  [m (/ (first ai) (first (first A)))
                   ai (map #(- %1 (* %2 m)) (rest ai) (rest (first A)))]
                  {:A (conj (AL :A) ai) :L (conj (AL :L) m)}))
              {:A [] :L []} (rest A)))
          (sub-ge [A k p]
            (when-first [a1 A]
              (let [Ap (pivot A 0)
                    Am (elimination-step (Ap :A))
                    li (reverse (Am :L))
                    ui (first (Ap :A))
                    pi (nth p (Ap :p))]
                (cons [li ui pi] ; TÄHÄN VIELÄ SUORA SIJOITUS MAP-RAKENTEESEEN!!!! JA RECUR- TAI REDUCE-RAKENNE
                      (sub-ge (Am :A) (inc k) (remove (partial = (Ap :p)) p))))))]
  (apply map vector (sub-ge A 0 (range (count A))))))


(def LUp (lu-factorization-sparse A))
(def LUp2 (lu-factorization-sparse A2))
(println "LUp = " LUp)
(println "LUp2 = " LUp2)

(defn lu-forward-substitution-sparse [L p b]
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
        b (range 1 (count b))))))

(def y (lu-forward-substitution-sparse (first LUp) (last LUp) [[8.3 -3.8 -2.3][5 6 2][1 8 1]]))
(println "y: " y)

(println "x: " (backward-substitution-sparse (second LUp) y))


(defn cholesky-decomposition [A]
  ; A must be symmetric posite-definite matrix.
  (reduce
    (fn [C i]
      (let [lower-units (reduce 
                          (fn [ci j]
                            (conj ci
	                            (/ 
	                              (- (nth (nth A i) j) (apply + (map * ci (nth C j))))
	                              (nth (nth C j) j 1))))
                          [] (range (count C)))
            diag-unit (Math/sqrt (- 
                                  (nth (nth A i) i) 
                                  (apply + (map #(* % %) lower-units))))]
        (conj C (conj lower-units diag-unit))))
    [] (range (count A))))
(cholesky-decomposition [[25 15 -5][15 18 0][-5 0 11]])
(cholesky-decomposition [[18 22 54 42][22 70 86 62][54 86 174 134][42 62 134 106]])
